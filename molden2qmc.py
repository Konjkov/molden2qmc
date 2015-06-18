#!/usr/bin/python
# -*- coding: utf-8 -*-

__version__ = '2.6.0'

"""
TODO:
1. implement Nwchem
2. implement QChem
"""

import os
from math import pi, sqrt, factorial, fabs
from itertools import combinations
from operator import itemgetter


def fact2(k):
    """
    Compute double factorial: k!! = 1*3*5*....k

    inspired by https://gist.github.com/fmeyer/289467
    """
    return reduce(float.__mul__, range(k, 0, -2), 1.0)


def smart_float(x):
    """
    Expect that x represent float in formats:
    '-12345.5678'
    '-123E+45'
    ' 123e-45'
    '-123D+45'
    ' 123d+45'
    """
    return float(x.replace('D', 'E').replace('d', 'e'))  # for MOLPRO only


def list_mul(list_a, list_b):
    """
    element-wise multiplication of two lists
    """
    return [a*b for a, b in zip(list_a, list_b)]


class Molden(object):
    """
    Data structures used in Molden Class:
    atom_list = [
        {'N': <atomic_number>,
         'X': <x-coordinate>,
         'Y': <y-coordinate>,
         'Z': <z-coordinate>,
         'pseudo': <True|False>,
         'SHELLS': [{'TYPE': <'s', 'p', 'd', ...>,
                     'DATA': [[exponent_primitive_1, contraction_coefficient_1],
                              [exponent_primitive_2, contraction_coefficient_2],
                              [exponent_primitive_3, contraction_coefficient_3],
                              ...
                             ]
                    },
                    shell2,
                    shell3,
                    ...
                   ]
        },
        atom2,
        atom3,
        ...
    ]
    mo_matrix = [
        {'SYMMETRY': <orbital symmetry>,
         'ENERGY': <orbital energy au>,
         'SPIN': <spin projection: alpha or beta>,
         'OCCUPATION': <0 or 1 or 2>,
         'MO': [{'TYPE': <'s', 'p', 'd', ...>,
                 'DATA': [mo_coefficient_1,
                          mo_coefficient_2,
                          ...
                         ] - list of (2l+1) coefficients
                },
                AO2,
                AO3,
                ...
               ]
        },
        mo_orbital2,
        mo_orbital3,
        ...
    ]
    """
    Ang2Bohr = 1/0.52917721  # 1 Bohr = 0.52917721 Angstrom
    ang_momentum_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'sp': 1}
    title = "Insert Your Title Here\n"

    def __init__(self, f, pseudoatoms="none"):
        """ Create instance that represent MOLDEN file data
        http://www.cmbi.ru.nl/molden/molden_format.html

        :param f: file descriptor of MOLDEN file
        :param pseudoatoms: list of pseudoatoms (none, all or white space separated numbers)
        """
        self.D_orb_conversion_required = True  # Conversion D orbitals Cartesian -> Spherical required
        self.F_orb_conversion_required = True  # Conversion F orbitals Cartesian -> Spherical required
        self.G_orb_conversion_required = True  # Conversion G orbitals Cartesian -> Spherical required
        self.pseudoatoms = pseudoatoms
        self.atom_list = []
        self.mo_matrix = []
        self.f = f
        self.molden_atoms()
        self.molden_pseudo()
        self.molden_gto()
        self.molden_mo()

    def molden_section(self, section_name):
        """
        :returns: content of named section
        """
        self.f.seek(0)
        line = self.f.readline()
        # CFOUR requires case insensitive comparison
        while line and not line.upper().startswith("[%s]" % section_name.upper()):
            line = self.f.readline()
        result = [line]
        line = self.f.readline()
        while line and not line.startswith('['):
            result.append(line)
            line = self.f.readline()
        return result

    def molden_atoms(self):
        """
        parse [Atoms] section.
        Format:
        [Atoms] (Angs|AU)
        element_name number atomic_number x y z
        """
        section = self.molden_section("Atoms")
        section_header = section[0]
        section_body = section[1:]

        for line in section_body:
            split_line = line.split()
            if section_header.split()[1] == 'Angs':
                atom = {'N': int(split_line[2]),  # atomic number
                        'X': float(split_line[3]) * self.Ang2Bohr,
                        'Y': float(split_line[4]) * self.Ang2Bohr,
                        'Z': float(split_line[5]) * self.Ang2Bohr}
            else:
                atom = {'N': int(split_line[2]),  # atomic number
                        'X': float(split_line[3]),
                        'Y': float(split_line[4]),
                        'Z': float(split_line[5])}
            self.atom_list.append(atom)

    def molden_pseudo(self):
        """
        set pseudopotential for atoms
        """

        if self.pseudoatoms == "none":
            for atom in self.atom_list:
                atom['pseudo'] = False
        elif self.pseudoatoms == "all":
            for atom in self.atom_list:
                atom['pseudo'] = True
        else:
            for i, atom in enumerate(self.atom_list):
                if str(i + 1) in self.pseudoatoms.split():
                    atom['pseudo'] = True
                else:
                    atom['pseudo'] = False

    def molden_gto(self):
        """
        parse [GTO] section.
        Format:
        [GTO]
        atom_sequence_number1 0
        shell_label number_of_primitives 1.00
        exponent_primitive_1 contraction_coefficient_1 (contraction_coefficient_1)
        ...
        empty line
        atom_sequence__number2 0
        shell_label number_of_primitives 1.00
        exponent_primitive_1 contraction_coefficient_1 (contraction_coefficient_1)
        ...
        empty line

        P.S.
        Lines with zero contraction coefficients are skipped.
        """
        section_body = self.molden_section("GTO")[1:]
        for line in section_body:
            split_line = line.split()
            if len(split_line) < 2:  # empty line
                pass
            elif len(split_line) == 2 and split_line[-1] == '0':  # new atom
                atom = self.atom_list[int(split_line[0])-1]
                atom['SHELLS'] = []
            elif len(split_line) == 3:   # new shell
                shell = {'TYPE': split_line[0], 'DATA': []}
                atom['SHELLS'].append(shell)
            else:
                if smart_float(split_line[1]) != 0.0:  # non zero contraction coefficient
                    shell['DATA'].append([smart_float(split_line[0]), smart_float(split_line[1])])

    def molden_spherical_cartesian(self):
        """
        Check that D, F, G orbitals required conversion from cartesian to spherical as described in documentation:
        Use the keyword [5D] on a separate line to specify the use of 'spherical' D and F functions
        (5 D and 7 F functions). The default is to use 'cartesian' D and F functions (6 D and 10 F functions).
        The enable the use of mixed spherical and cartesian function, the following keywords where added
        ([5D10F], [7F] (6D en 7F), [5D7F], (same as[5D], for reasons of backwards compatibility).
        Since molden 4.4 G-functions are also supported, default is cartesian G functions.
        Use [9G] to specify spherical G functions.

        Conversion required by default.
        """
        self.f.seek(0)
        for line in self.f:
            line = line.upper()  # PSI4 requires case insensitive comparison
            if line.startswith("[5D]") or line.startswith("[5D7F]"):
                self.D_orb_conversion_required = False
                self.F_orb_conversion_required = False
            if line.startswith("[5D10F]"):
                self.D_orb_conversion_required = False
                self.F_orb_conversion_required = True
            if line.startswith("[7F]"):
                self.D_orb_conversion_required = True
                self.F_orb_conversion_required = False
            if line.startswith("[9G]"):
                self.G_orb_conversion_required = False

    def molden_mo(self):
        """
        parse [MO] section.
        Format:
        [MO]
        Sym= symmetry_label_1
        Ene= mo_energy_1
        Spin= (Alpha|Beta)
        Occup= mo_occupation_number_1
        ao_number_1 mo_coefficient_1
        ...
        ao_number_n mo_coefficient_n
        ....
        Sym= symmetry_label_N
        Ene= mo_energy_N
        Spin= (Alpha|Beta)
        Occup= mo_occupation_number_N
        ao_number_1 mo_coefficient_1
        ...
        ao_number_n mo_coefficient_n
        """
        self.molden_spherical_cartesian()
        # (Number of basis functions) blocks of (Number of basis functions) lines each
        mo_length_map = {'s': 1,
                         'p': 3,
                         'd': 6 if self.D_orb_conversion_required else 5,
                         'f': 10 if self.F_orb_conversion_required else 7,
                         'g': 15 if self.G_orb_conversion_required else 9}

        nbasisfunctions = 0
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                    nbasisfunctions += mo_length_map[shell['TYPE']]

        section_body = self.molden_section("MO")[1:]
        for line in section_body:
            if line.strip().startswith('Sym='):
                # Dalton remove lines with zero coefficients, restoring them
                mo_orbital_block = {'raw_data': [0.0] * nbasisfunctions}
                self.mo_matrix.append(mo_orbital_block)
                mo_orbital_block['SYMMETRY'] = line.split('=')[1].strip()
            elif line.strip().startswith('Ene='):
                mo_orbital_block['ENERGY'] = float(line.split('=')[1])
            elif line.strip().startswith('Spin='):
                # ORCA don't put a space between SPIN= and Beta.
                mo_orbital_block['SPIN'] = line.split('=')[1].strip()
            elif line.strip().startswith('Occup='):
                mo_orbital_block['OCCUPATION'] = float(line.split('=')[1])
            else:
                split_line = line.split()
                mo_orbital_block['raw_data'][int(split_line[0])-1] = float(split_line[1])

        for mo_orbital_block in self.mo_matrix:
            offset = 0
            mo_orbital_block['MO'] = []
            for atom in self.atom_list:
                for shell in atom['SHELLS']:
                    shell_length = mo_length_map[shell['TYPE']]
                    ao = {'TYPE': shell['TYPE'],
                          'DATA': mo_orbital_block['raw_data'][offset:offset+shell_length]}
                    mo_orbital_block['MO'].append(ao)
                    offset += shell_length

    def spin_unrestricted(self):
        """
        :return: True when wfn is unrestricted or False elsewhere.
        """
        return reduce(bool.__or__, [(mo_orbital_block['SPIN'] == 'Beta')
                                    for mo_orbital_block in self.mo_matrix], False)

    def charge(self, atom):
        """
        :returns: valence charge of atom.
        """
        if atom['pseudo']:
            if atom['N'] <= 2:     # H-He
                return atom['N']
            elif atom['N'] <= 10:  # Li-Ne
                return atom['N'] - 2
            elif atom['N'] <= 18:  # Na-Ar
                return atom['N'] - 10
            elif atom['N'] <= 26:  # K-Zn
                return atom['N'] - 18
            elif atom['N'] <= 36:  # Ga-Kr
                return atom['N'] - 28
            elif atom['N'] <= 48:  # Rb-Cd
                return atom['N'] - 36
            elif atom['N'] <= 54:  # In-Xe
                return atom['N'] - 46
            elif atom['N'] <= 71:  # Cs-Lu
                return atom['N'] - 54
            elif atom['N'] <= 80:  # Hf-Hg
                return atom['N'] - 68
            else:
                raise NotImplementedError("AREP Trail & Needs PP didn't support elements after Hg")
        else:
            return atom['N']

    def natom(self):
        """
        :returns: total number of atoms
        """
        return len(self.atom_list)

    def nelec(self):
        """
        :returns: total number of electrons
        """
        return sum(orbital['OCCUPATION'] for orbital in self.mo_matrix)

    def nshell(self):
        """
        :returns: total number of shells
        """
        return sum(len(atom['SHELLS']) for atom in self.atom_list)

    def nbasis_functions(self):
        """
        :returns: total number of basis functions converted to spherical
        """
        result = 0
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                    result += 2 * self.ang_momentum_map[shell['TYPE']] + 1
        return result

    def nprimitives(self):
        """
        :returns: total number of primitives
        """
        result = 0
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                result += len(shell['DATA'])
        return result

    def highest_ang_mo(self):
        """
        (s/p/d/f... 1/2/3/4...)
        :returns: highest angular momentum
        """
        result = 0
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                result = max(result, self.ang_momentum_map[shell['TYPE']] + 1)
        return result

    def distance(self, atom1, atom2):
        """
        :returns: distance between atoms
        """
        return sqrt((atom1['X'] - atom2['X'])**2 +
                    (atom1['Y'] - atom2['Y'])**2 +
                    (atom1['Z'] - atom2['Z'])**2)

    def nuclear_repulsion(self):
        """
        :returns: n-n repulsion energy
        """
        return sum(self.charge(atom1) * self.charge(atom2)/self.distance(atom1, atom2)
                   for atom1, atom2 in combinations(self.atom_list, 2))

    def gwfn(self):
        """
        write out gwfn.data file
        """
        gwfn = open('gwfn.data', 'w')
        gwfn.write(self.gwfn_title())
        gwfn.write(self.gwfn_basic_info())
        gwfn.write(self.gwfn_geometry())
        gwfn.write(self.gwfn_basis_set())
        gwfn.write(self.gwfn_multideterminant_information())
        gwfn.write(self.gwfn_orbital_coefficients())

    def gwfn_title(self):
        """
        :returns: TITLE section of gwfn.data file
        """
        return ("TITLE\n"
                "%s") % self.title

    def gwfn_basic_info(self):
        """
        :returns: BASIC_INFO section of gwfn.data file
        """
        return ("BASIC_INFO\n"
                "----------\n"
                "Generated by:\n"
                "    molden2qmc version: %s\n"
                "Method:\n"
                "\n"
                "DFT Functional:\n"
                "\n"
                "Periodicity:\n"
                "          0\n"
                "Spin unrestricted:\n"
                "          %s\n"
                "nuclear-nuclear repulsion energy (au/atom):\n"
                "          %s\n"
                "Number of electrons per primitive cell:\n"
                "          %s\n"
                "\n") % (__version__,
                         '.true.' if self.spin_unrestricted() else '.false.',
                         self.nuclear_repulsion()/self.natom(),
                         int(round(self.nelec())))

    def gwfn_geometry(self):
        """
        :returns: GEOMETRY section of gwfn.data file
        """
        result = ("GEOMETRY\n"
                  "--------\n"
                  "Number of atoms:\n"
                  "          %s\n"
                  "Atomic positions (au):\n") % self.natom()

        for atom in self.atom_list:
            result += "% .13E% .13E% .13E\n" % (atom['X'], atom['Y'], atom['Z'])

        result += "Atomic numbers for each atom:"
        for num, atom in enumerate(self.atom_list):
            if num % 8 == 0:
                result += "\n"
            if atom['pseudo']:
                result += "%10d" % (atom['N'] + 200)
            else:
                result += "%10d" % atom['N']

        result += '\nValence charges for each atom:'
        for num, atom in enumerate(self.atom_list):
            if num % 4 == 0:
                result += "\n"
            result += " %1.13E" % self.charge(atom)

        return result + "\n\n"

    def gwfn_basis_set(self):
        """
        :returns: BASIS SET section of gwfn.data file
        """
        result = ("BASIS SET\n"
                  "---------\n"
                  "Number of Gaussian centres\n"
                  "         %s\n"
                  "Number of shells per primitive cell\n"
                  "         %s\n"
                  "Number of basis functions (\'AO\') per primitive cell\n"
                  "         %s\n"
                  "Number of Gaussian primitives per primitive cell\n"
                  "         %s\n"
                  "Highest shell angular momentum (s/p/d/f... 1/2/3/4...)\n"
                  "         %s\n"
                  "Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)") % (self.natom(),
                                                                          self.nshell(),
                                                                          self.nbasis_functions(),
                                                                          self.nprimitives(),
                                                                          self.highest_ang_mo())
        num = 0
        ang_type_map = {'s': 1, 'sp': 2, 'p': 3, 'd': 4, 'f': 5, 'g': 6}
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                if num % 8 == 0:
                    result += "\n"
                result += "%10d" % ang_type_map[shell['TYPE']]
                num += 1

        num = 0
        result += "\nNumber of primitive Gaussians in each shell"
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                if num % 8 == 0:
                    result += "\n"
                result += "%10d" % len(shell['DATA'])
                num += 1

        sequence_number = 1
        dummy_atom = [{'SHELLS': ()}]  # hack
        result += "\nSequence number of first shell on each centre"
        for num, atom in enumerate(dummy_atom + self.atom_list):
            sequence_number += len(atom['SHELLS'])
            if num % 8 == 0:
                result += "\n"
            result += "%10d" % sequence_number

        num = 0
        result += "\nExponents of Gaussian primitives"
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                for primitive in shell['DATA']:
                    if num % 4 == 0:
                        result += "\n"
                    result += "% .13E" % primitive[0]
                    num += 1

        num = 0
        result += "\nNormalized contraction coefficients"
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                for primitive in shell['DATA']:
                    if num % 4 == 0:
                        result += "\n"
                    result += "% .13E" % primitive[1]
                    num += 1

        result += "\nPosition of each shell (au)"
        for atom in self.atom_list:
            for _ in atom['SHELLS']:  # throwaway variable
                result += "\n% .13E% .13E% .13E" % (atom['X'], atom['Y'], atom['Z'])

        return result + "\n\n"

    def gwfn_multideterminant_information(self):
        """
        :returns: MULTIDETERMINANT INFORMATION section of gwfn.data file
        """
        return ("MULTIDETERMINANT INFORMATION\n"
                "----------------------------\n"
                "GS\n\n")

    def gwfn_orbital_coefficients(self):
        """
        :returns: ORBITAL COEFFICIENTS section of gwfn.data file, sorted by ('SPIN', 'ENERGY')
        """
        result = ("ORBITAL COEFFICIENTS\n"
                  "------------------------")

        # (Number of basis functions) ** 2 coefficients
        num = 0
        for orbital in sorted(self.mo_matrix, key=itemgetter('SPIN', 'ENERGY')):
            for ao in orbital['MO']:
                for coefficient in ao['DATA']:
                    if num % 4 == 0:
                        result += "\n"
                    result += "% .13E" % coefficient
                    num += 1

        return result + "\n\n"


class DefaultConverter(Molden):
    """
    Default converter expect that contraction coefficients is
    'Published' in the EMSL Basis Set Library, MO-coefficients is
    in a spherical format and only m-dependent normalisation required.
    """

    def __init__(self, f, pseudoatoms="none"):
        super(DefaultConverter, self).__init__(f, pseudoatoms)
        self.atom_list_converter()
        self.mo_matrix_converter()

    def whole_contraction_factor(self, primitives, l):
        """
        :param primitives: ((a_1, d_1), (a_2, d_2), ...)
        :param l: angular quantum number

        The normalization constant N_cont for the whole contraction (in which the
        contraction coefficient for the ith primitive is d_i, and the exponent is a_i)
        is given by:
                                                1
        N_cont = ---------------------------------------------------------
                                         2 * root (a_i * a_j ).
                 sqrt [ sum_ij d_i d_j ( --------------------- )^(l+3/2) ]
                                               a_i + a_j
        """
        s = 0
        for p1 in primitives:
            for p2 in primitives:
                s += p1[1] * p2[1] * (2 * sqrt(p1[0] * p2[0])/(p1[0] + p2[0]))**(l + 1.5)
        return 1/sqrt(s)

    def m_independent_factor(self, a, l):
        """
        :param a: alpha
        :param l: angular quantum number

        The m-independent factors for the different shells:
                 root[2^(l+3/2) * alpha^(l+3/2)]          2^l
        N_prim = ------------------------------- * root --------
                        pi^(3/4)                        (2l-1)!!
        """
        return sqrt(2**(l + 1.5) * a**(l + 1.5))/pi**0.75 * sqrt(2**l/fact2(2*l-1))

    def m_dependent_factor(self, l, m):
        """
        :param l: angular quantum number
        :param m: magnetic quantum number
        The m-dependent factors for the different shells:

               (2 - delta_m,0) * (l - |m|)!
        root  ------------------------------
                     (l + |m|)!

        read examples/generic/gauss_dfg/README for details
        """
        if l < 2 or m == 0:
            return 1.0
        else:
            return sqrt(2.0 * factorial(l - fabs(m))/factorial(l + fabs(m)))

    def atom_list_converter(self):
        """
        Default contraction coefficients is 'Published' in the EMSL Basis Set Library.
        """
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                l = self.ang_momentum_map[shell['TYPE']]
                w = self.whole_contraction_factor(shell['DATA'], l)
                for primitive in shell['DATA']:
                    primitive[1] *= w * self.m_independent_factor(primitive[0], l)

    def d_to_spherical(self, cartesian):
        """
        Convert cartesian representation of d-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of D functions is expected:
            5D: D 0, D+1, D-1, D+2, D-2
            6D: xx, yy, zz, xy, xz, yz
        """
        xx, yy, zz, xy, xz, yz = cartesian

        r2 = xx + yy + zz

        zero = (3.0 * zz - r2) / 2.0
        plus_1 = sqrt(3) * xz
        minus_1 = sqrt(3) * yz
        plus_2 = sqrt(3) * (xx - yy) / 2.0
        minus_2 = sqrt(3) * xy
        return zero, plus_1, minus_1, plus_2, minus_2

    def f_to_spherical(self, cartesian):
        """
        Convert cartesian representation of f-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of F functions is expected:
            7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
            10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        """
        xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz = cartesian

        xr2 = xxx + xyy + xzz
        yr2 = xxy + yyy + yzz
        zr2 = xxz + yyz + zzz

        zero = (5.0 * zzz - 3.0 * zr2) / 2.0
        plus_1 = sqrt(6) * (5.0 * xzz - xr2) / 4.0
        minus_1 = sqrt(6) * (5.0 * yzz - yr2) / 4.0
        plus_2 = sqrt(15) * (xxz - yyz) / 2.0
        minus_2 = sqrt(15) * xyz
        plus_3 = sqrt(10) * (xxx - 3.0 * xyy) / 4.0
        minus_3 = sqrt(10) * (3.0 * xxy - yyy) / 4.0
        return zero, plus_1, minus_1, plus_2, minus_2, plus_3, minus_3

    def g_to_spherical(self, cartesian):
        """
        Convert cartesian representation of g-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of G functions is expected:
            9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
            15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
                 xxyy xxzz yyzz xxyz yyxz zzxy
        """
        xxxx, yyyy, zzzz, xxxy, xxxz, yyyx, yyyz, zzzx, zzzy, xxyy, xxzz, yyzz, xxyz, yyxz, zzxy = cartesian

        xyr2 = xxxy + yyyx + zzxy
        xzr2 = xxxz + yyxz + zzzx
        yzr2 = xxyz + yyyz + zzzy
        x2r2 = xxxx + xxyy + xxzz
        y2r2 = xxyy + yyyy + yyzz
        z2r2 = xxzz + yyzz + zzzz
        r4 = x2r2 + y2r2 + z2r2

        zero = (35.0 * zzzz - 30.0 * z2r2 + 3.0 * r4) / 8.0
        plus_1 = sqrt(10) * (7.0 * zzzx - 3.0 * xzr2) / 4.0
        minus_1 = sqrt(10) * (7.0 * zzzy - 3.0 * yzr2) / 4.0
        plus_2 = sqrt(5) * (7.0 * (xxzz - yyzz) - (x2r2 - y2r2)) / 4.0
        minus_2 = sqrt(5) * (7.0 * zzxy - xyr2) / 2.0
        plus_3 = sqrt(70) * (xxxz - 3.0 * yyxz) / 4.0
        minus_3 = sqrt(70) * (3.0 * xxyz - yyyz) / 4.0
        plus_4 = sqrt(35) * (xxxx - 6.0 * xxyy + yyyy) / 8.0
        minus_4 = sqrt(35) * (xxxy - yyyx) / 2.0
        return zero, plus_1, minus_1, plus_2, minus_2, plus_3, minus_3, plus_4, minus_4

    def d_normalize(self, coefficient):
        """
        The following order of D functions is expected:
            5D: D 0, D+1, D-1, D+2, D-2

        P.S.
        One historical CASINO inconsistency which may be easily overlooked:
        Constant numerical factors in the real solid harmonics e.g. the '3' in the 3xy
        d function, or '15' in the (15x^3-45^xy2) f function, may be premultiplied into
        the orbital coefficients so that CASINO doesn't have to e.g. multiply by 3
        every time it evaluates that particular d function. In practice the CASINO
        orbital evaluators do this only for d functions, but *not for f and g* (this
        may or may not be changed in the future if it can be done in a.
        backwards-consistent way)
        """
        premultiplied_factor = (0.5, 3.0, 3.0, 3.0, 6.0)
        return (coefficient[0] * self.m_dependent_factor(2,  0) * premultiplied_factor[0],
                coefficient[1] * self.m_dependent_factor(2,  1) * premultiplied_factor[1],
                coefficient[2] * self.m_dependent_factor(2, -1) * premultiplied_factor[2],
                coefficient[3] * self.m_dependent_factor(2,  2) * premultiplied_factor[3],
                coefficient[4] * self.m_dependent_factor(2, -2) * premultiplied_factor[4])

    def f_normalize(self, coefficient):
        """
        The following order of F functions is expected:
            7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
        """
        return (coefficient[0] * self.m_dependent_factor(3,  0),
                coefficient[1] * self.m_dependent_factor(3,  1),
                coefficient[2] * self.m_dependent_factor(3, -1),
                coefficient[3] * self.m_dependent_factor(3,  2),
                coefficient[4] * self.m_dependent_factor(3, -2),
                coefficient[5] * self.m_dependent_factor(3,  3),
                coefficient[6] * self.m_dependent_factor(3, -3))

    def g_normalize(self, coefficient):
        """
        The following order of G functions is expected:
            9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
        """
        return (coefficient[0] * self.m_dependent_factor(4,  0),
                coefficient[1] * self.m_dependent_factor(4,  1),
                coefficient[2] * self.m_dependent_factor(4, -1),
                coefficient[3] * self.m_dependent_factor(4,  2),
                coefficient[4] * self.m_dependent_factor(4, -2),
                coefficient[5] * self.m_dependent_factor(4,  3),
                coefficient[6] * self.m_dependent_factor(4, -3),
                coefficient[7] * self.m_dependent_factor(4,  4),
                coefficient[8] * self.m_dependent_factor(4, -4))

    def mo_matrix_converter(self):
        """
        Only mo_coefficients of d, f, g must be converted by default.
        """
        for orbital in self.mo_matrix:
            for ao in orbital['MO']:
                if ao['TYPE'] == 'd':
                    if self.D_orb_conversion_required:
                        ao['DATA'] = self.d_to_spherical(ao['DATA'])
                    ao['DATA'] = self.d_normalize(ao['DATA'])
                elif ao['TYPE'] == 'f':
                    if self.F_orb_conversion_required:
                        ao['DATA'] = self.f_to_spherical(ao['DATA'])
                    ao['DATA'] = self.f_normalize(ao['DATA'])
                elif ao['TYPE'] == 'g':
                    if self.G_orb_conversion_required:
                        ao['DATA'] = self.g_to_spherical(ao['DATA'])
                    ao['DATA'] = self.g_normalize(ao['DATA'])


class Turbomole(DefaultConverter):
    """
    Turbomole 6.6
    """
    title = "generated from Turbomole output data.\n"

    def d_to_spherical(self, cartesian):
        """
        Convert cartesian representation of d-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of D functions is expected:
            5D: D 0, D+1, D-1, D+2, D-2
            6D: xx, yy, zz, xy, xz, yz
        """
        norm = [2.0/sqrt(3)] * 3 + [1.0] * 3
        return super(Turbomole, self).d_to_spherical(list_mul(norm, cartesian))

    def f_to_spherical(self, cartesian):
        """
        Convert cartesian representation of f-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of F functions is expected:
            7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
            10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        """
        norm = [6.0/sqrt(15)] * 3 + [2.0/sqrt(3)] * 6 + [1.0]
        return super(Turbomole, self).f_to_spherical(list_mul(norm, cartesian))

    def g_to_spherical(self, cartesian):
        """
        Convert cartesian representation of g-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of G functions is expected:
            9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
            15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
                 xxyy xxzz yyzz xxyz yyxz zzxy
        """
        norm = [24.0/sqrt(105)] * 3 + [6.0/sqrt(15)] * 6 + [4.0/3.0] * 3 + [2.0/sqrt(3)] * 3
        return super(Turbomole, self).g_to_spherical(list_mul(norm, cartesian))


class CFour(DefaultConverter):
    """
    CFour 2.0 beta
    In order to get the correct MOLDEN file recompile Cfour with modified reorderedf.f
    https://github.com/Konjkov/molden2qmc/blob/master/test/CFOUR/reorderdf.f
    """
    title = "generated from CFour output data.\n"

    def nelec(self):
        """
        in CFOUR occupation number sometimes takes value from list (0, 1)
        """
        if self.spin_unrestricted():
            return super(CFour, self).nelec()
        else:
            return 2 * super(CFour, self).nelec()

    def d_to_spherical(self, cartesian):
        """
        Convert cartesian representation of d-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of D functions is expected:
            5D: D 0, D+1, D-1, D+2, D-2
            6D: xx, yy, zz, xy, xz, yz
        """
        norm = [2.0/sqrt(3)] * 3 + [1.0/sqrt(3)] * 3
        return super(CFour, self).d_to_spherical(list_mul(norm, cartesian))

    def f_to_spherical(self, cartesian):
        """
        Convert cartesian representation of f-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of F functions is expected:
            7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
            10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        """
        norm = [6.0/sqrt(15)] * 3 + [2.0/sqrt(15)] * 6 + [1.0/sqrt(15)]
        return super(CFour, self).f_to_spherical(list_mul(norm, cartesian))

    def g_to_spherical(self, cartesian):
        """
        Convert cartesian representation of g-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of G functions is expected:
            9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
            15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
                 xxyy xxzz yyzz xxyz yyxz zzxy
        """
        norm = [24.0/sqrt(105)] * 3 + [6.0/sqrt(105)] * 6 + [4.0/sqrt(105)] * 3 + [2.0/sqrt(105)] * 3
        return super(CFour, self).g_to_spherical(list_mul(norm, cartesian))


class Orca(DefaultConverter):
    """
    ORCA 3.X
    """
    title = "generated from Orca output data.\n"

    def d_to_spherical(self, cartesian):
        """
        In ORCA spherical MOLDEN input expected.
        """
        raise RuntimeWarning("Cartesian input is not expected. Check that"
                             " the quantum chemistry code is selected correctly")

    g_to_spherical = f_to_spherical = d_to_spherical

    def atom_list_converter(self):
        """
        in ORCA 's', 'p' orbitals don't require normalization.
        'g' orbitals need to be additionally scaled up by a factor of sqrt(3).
        https://orcaforum.cec.mpg.de/viewtopic.php?f=8&t=1484
        """
        for atom in self.atom_list:
            for shell in atom['SHELLS']:
                l = self.ang_momentum_map[shell['TYPE']]
                for primitive in shell['DATA']:
                    primitive[1] /= sqrt(fact2(2*l-1))
                    if l == 4:
                        primitive[1] *= sqrt(3)

    def f_normalize(self, coefficient):
        """
        ORCA use slightly different sign conventions:
            F(+3)_ORCA = - F(+3)_MOLDEN
            F(-3)_ORCA = - F(-3)_MOLDEN
        """
        coefficient[5] *= -1
        coefficient[6] *= -1
        return super(Orca, self).f_normalize(coefficient)

    def g_normalize(self, coefficient):
        """
        ORCA use slightly different sign conventions:
            G(+3)_ORCA = - G(+3)_MOLDEN
            G(-3)_ORCA = - G(-3)_MOLDEN
            G(+4)_ORCA = - G(+4)_MOLDEN
            G(-4)_ORCA = - G(-4)_MOLDEN
        """
        coefficient[5] *= -1
        coefficient[6] *= -1
        coefficient[7] *= -1
        coefficient[8] *= -1
        return super(Orca, self).g_normalize(coefficient)


class PSI4(DefaultConverter):
    """
    PSI4 https://github.com/psi4/psi4public
    """
    title = "generated from PSI4 output data.\n"

    def d_to_spherical(self, cartesian):
        """
        In PSI4 spherical MOLDEN input expected.
        """
        raise RuntimeWarning("Cartesian input is not expected. Check that"
                             " the quantum chemistry code is selected correctly")

    g_to_spherical = f_to_spherical = d_to_spherical

    def atom_list_converter(self):
        """
        PSI4 correctly normalized contraction coefficients
        in the MOLDEN file, so we don't need to modify them.
        """


class Dalton(DefaultConverter):
    """
    DALTON 2013
    """
    title = "generated from Dalton output data.\n"

    def d_to_spherical(self, cartesian):
        """
        In DALTON spherical MOLDEN input expected.
        """
        raise RuntimeWarning("Cartesian input is not expected. Check that"
                             " the quantum chemistry code is selected correctly")

    g_to_spherical = f_to_spherical = d_to_spherical


class Molpro(DefaultConverter):
    """
    Molpro
    """
    title = "generated from Molpro output data.\n"

    def d_to_spherical(self, cartesian):
        """
        Convert cartesian representation of d-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of D functions is expected:
            5D: D 0, D+1, D-1, D+2, D-2
            6D: xx, yy, zz, xy, xz, yz
        """
        norm = [2.0/sqrt(3)/sqrt(3)] * 3 + [1.0/sqrt(3)] * 3
        return super(Molpro, self).d_to_spherical(list_mul(norm, cartesian))

    def f_to_spherical(self, cartesian):
        """
        Convert cartesian representation of f-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of F functions is expected:
            7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
            10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        """
        norm = [6.0/sqrt(15)/sqrt(15)] * 3 + [2.0/sqrt(3)/sqrt(15)] * 6 + [1.0/sqrt(15)]
        return super(Molpro, self).f_to_spherical(list_mul(norm, cartesian))

    def g_to_spherical(self, cartesian):
        """
        Convert cartesian representation of g-orbital to spherical
        http://theochem.github.io/horton/tut_gaussian_basis.html
        The following order of G functions is expected:
            9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
            15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
                 xxyy xxzz yyzz xxyz yyxz zzxy
        """
        norm = [24.0/sqrt(105)/sqrt(105)] * 3 + [6.0/sqrt(15)/sqrt(105)] * 6 +\
               [4.0/3.0/sqrt(105)] * 3 + [2.0/sqrt(3)/sqrt(105)] * 3

        return super(Molpro, self).g_to_spherical(list_mul(norm, cartesian))


def main():
    print ("Hello, you are converting a MOLDEN file to a CASINO gwfn.data file.\n")

    while True:
        input_file_name = raw_input("Enter the name of your MOLDEN file: ")
        if os.path.exists(str(input_file_name)):
            break
        print ("File not found...")

    input_file = open(str(input_file_name), "r")

    code = raw_input("\nEnter the NUMBER corresponding to the quantum chemistry "
                     "code used to produce this MOLDEN file:\n"
                     "0 -- TURBOMOLE\n"
                     "1 -- PSI4\n"
                     "2 -- CFOUR 2.0beta\n"
                     "3 -- ORCA 3.X\n"
                     "4 -- DALTON2013\n"
                     "5 -- MOLPRO\n")
    while not code.isdigit() or int(code) > 5:
        code = raw_input('Sorry,  try again.')

    code = int(code)
    print "You have entered NUMBER ", code, "\n"

    pseudoatoms = raw_input("This script did not detect if a pseudopotential was used.\n"
                            "Please eneter the list of atoms for thouse pseudopotential was used:\n\n"
                            "none = a pseudopotential was not used for any atoms in this calculation.\n"
                            "all = a pseudopotential was used for all atoms in this calculation.\n"
                            "white space separated numbers = number of pseudoatoms (started from 1).\n\n").lower()

    print "You have entered ", pseudoatoms, "\n"

    if code == 0:
        Turbomole(input_file, pseudoatoms).gwfn()
    elif code == 1:
        PSI4(input_file, pseudoatoms).gwfn()
    elif code == 2:
        CFour(input_file, pseudoatoms).gwfn()
    elif code == 3:
        Orca(input_file, pseudoatoms).gwfn()
    elif code == 4:
        Dalton(input_file, pseudoatoms).gwfn()
        print ("In Dalton's MOLDEN file all occupation numbers of HF and DFT MOs are zero values by mistake\n"
               "So you should correct 'Number of electrons per primitive cell' in gwfn.data file by hand.")
    elif code == 5:
        Molpro(input_file, pseudoatoms).gwfn()

if __name__ == "__main__":
    main()
