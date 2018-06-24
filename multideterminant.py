#!/usr/bin/env python3
import re
import sys
import argparse
from functools import reduce


class SectionNotFound(Exception):
    """Section not found in MOLDEN file."""

    def __init__(self, section_name):
        self.section_name = section_name

    def __str__(self):
        return repr(self.section_name)


class QChemSectionNotFound(Exception):
    """Section not found in QChem Output file."""


class ORCASectionNotFound(Exception):
    """Section not found in ORCA Output file."""


class Default:
    """All possible combinations."""

    title = "All possible combinations of DETs."

    def __init__(self):
        """Initialise multi-determinant support."""
        self.max_orbital = 18
        self.alpha_occ = (6, 7)
        self.beta_occ = (6,)
        self.alpha_virt = range(6, self.max_orbital + 1)
        self.beta_virt = range(6, self.max_orbital + 1)
        self.orbitals = (
                len(self.alpha_occ) *
                len(self.beta_occ) *
                len(self.alpha_virt) *
                len(self.beta_virt)
        ) + 1

    def correlation(self):
        """
        :returns: MDET section of correlation.data file
        """
        with open('correlation.data', 'w') as output_file:
            print('START MDET', file=output_file)
            print('Title', file=output_file)
            print(' multideterminant WFN %s\n' % self.title, file=output_file)

            print('MD', file=output_file)
            print('  %i' % self.orbitals, file=output_file)
            print('  1.0  1  0', file=output_file)
            for i in range(1, self.orbitals):
                print('  0.0  %i  1' % (i+1), file=output_file)

            i = 1
            for alpha_from in self.alpha_occ:
                for beta_from in self.beta_occ:
                    for alpha_to in self.alpha_virt:
                        for beta_to in self.beta_virt:
                            i += 1
                            print('  DET %i 1 PR %i 1 %i 1' % (i, alpha_from, alpha_to), file=output_file)
                            print('  DET %i 2 PR %i 1 %i 1' % (i, beta_from, beta_to), file=output_file)
            print('END MDET', file=output_file)


class QChem:
    """
    QChem 4.4
    """
    title = "generated from QChem output data."

    def __init__(self, output_file, order=0, tol=0.01):
        """Initialise multi-determinant support."""
        self.output_file = output_file
        self.occupied = {'alpha': 0, 'beta': 0}
        self.active = {'alpha': 0, 'beta': 0}
        self.determinants = []
        self.parse_output()
        self.truncate(order, tol)

    def parse_output(self):
        """Retrieve from QChem output:
            T2-amplitudes information form VOD output in following format:
        -0.1500      5( Ag  ) A,   5( Ag  ) B  ->   1( B2u ) A,   1( B2u ) B  (VV)
        :return: list of spin-determinants, number of internal orbitals
        """
        electron_regexp = re.compile(
            'There are\s+(?P<alpha>\d+) alpha '
            'and\s+(?P<beta>\d+) beta electrons'
        )
        amplitude_regexp = re.compile(
            '(?P<weight>[-+]?\d+\.\d+)\s+'
            '(?P<first_from>\d+)\(.{5}\) '
            '(?P<first_from_spin>[AB]),\s+'
            '(?P<second_from>\d+)\(.{5}\) '
            '(?P<second_from_spin>[AB])  ->\s+'
            '(?P<first_to>\d+)\(.{5}\) '
            '(?P<first_to_spin>[AB]),\s+'
            '(?P<second_to>\d+)\(.{5}\) '
            '(?P<second_to_spin>[AB])'
        )
        spin_map = {'A': 1, 'B': 2}
        with open(self.output_file, "r") as qchem_output:
            line = qchem_output.readline()
            while line and not line.startswith('       Value      i             j           ->   a             b'):
                m = re.search(electron_regexp, line)
                if m:
                    self.occupied['alpha'] = int(m.group('alpha'))
                    self.occupied['beta'] = int(m.group('beta'))
                line = qchem_output.readline()
            if not line:
                raise QChemSectionNotFound('T2-amplitudes')
            line = qchem_output.readline()
            virtual_map = {'A': self.occupied['alpha'], 'B': self.occupied['beta']}
            while line and not line == '\n':
                m = re.search(amplitude_regexp, line)
                self.determinants.append({
                    'weight': float(m.group('weight')),
                    'promotions': (
                        {'from': int(m.group('first_from')) + 1,
                         'spin': spin_map[m.group('first_from_spin')],
                         'to': int(m.group('first_to')) + 1 + virtual_map[m.group('first_from_spin')]}
                        ,
                        {'from': int(m.group('second_from')) + 1,
                         'spin': spin_map[m.group('second_from_spin')],
                         'to': int(m.group('second_to')) + 1 + virtual_map[m.group('second_from_spin')]}
                    )
                })
                line = qchem_output.readline()

    def truncate(self, order=None, tol=0.01):
        """Leave only determinants with active space orbital
        number not greater then order."""
        determinants = []
        limit = order + self.occupied['alpha']
        for det in self.determinants:
            if order and (det['promotions'][0]['to'] > limit or det['promotions'][1]['to'] > limit):
                continue
            if abs(det['weight']) < tol:
                break
            determinants.append(det)
        self.determinants = determinants

    def correlation(self):
        """
        :returns: MDET section of correlation.data file
        """
        if len(self.determinants) > 0:
            with open('correlation.data', 'w') as output_file:
                print('START MDET', file=output_file)
                print('Title', file=output_file)
                print(' multideterminant WFN %s\n' % self.title, file=output_file)

                print('MD', file=output_file)
                print('  %i' % (len(self.determinants) + 1), file=output_file)
                # first determinant
                print(' %9.6f    1    0' % 1.0, file=output_file)
                for i, det in enumerate(self.determinants):
                    print(' %9.6f    %i    1' % (det['weight'], i+2), file=output_file)

                for i, det in enumerate(self.determinants):
                    for p in det['promotions']:
                        # starting from 2-nd determinant
                        print('  DET %i %i PR %i 1 %i 1' % (i + 2, p['spin'], p['from'], p['to']), file=output_file)
                print('END MDET', file=output_file)


class Orca:
    """
    ORCA 4.X
    """
    title = "generated from Orca output data"

    def __init__(self, orca_input_path):
        """Initialise multi-determinant support."""
        self.internal = 0  # CASSCF internal orbitals
        self.active = 0  # CASSCF active orbitals
        self.spin_determinants = []  # CASSCF spin-determinants
        self.orca_input_path = orca_input_path
        self.parse_output()
        # self.check_monotonic()

    def orca_section(self, section_name):
        """
        :returns: content of named section
        """
        with open(self.orca_input_path, "r") as orca_input:
            orca_input.seek(0)
            line = orca_input.readline()
            while line and not line.startswith(section_name):
                line = orca_input.readline()
            if not line:
                raise SectionNotFound(section_name)
            result = [line]
            line = orca_input.readline()
            line = orca_input.readline()
            while line and not line.startswith('-------'):
                result.append(line)
                line = orca_input.readline()
            return result

    def parse_output(self):
        """Retrive from ORCA output:
            Spin-Determinant information
            number of active & internal orbitals in CASSCF

        because of this issue https://orcaforum.cec.mpg.de/viewtopic.php?f=8&t=3212
        ORCA input should be:
        ! CASSCF cc-pVQZ

        %casscf
          nel  3
          norb 4
          PrintWF det
        end

        * xyzfile 0 2 ../../mol.xyz

        $new_job
        %casscf
          nel  3
          norb 4
          PrintWF csf
        end

        * xyzfile 0 2 ../../mol.xyz

        Keep your attention on not using symmetry because of this issue
        https://orcaforum.cec.mpg.de/viewtopic.php?f=11&t=3254

        :return: list of spin-determinants, number of internal orbitals
        """

        determinants = {}

        def ternary2decimal(ternary):
            """Convert ternary string to decimal.
            '001122' = 1 * 2 + 3 * 2 + 9 * 1 + 27 * 1 = 44
            '21102' -> 200
            """
            return reduce(lambda x, y: x * 3 + int(y), ternary, 0)

        with open(self.orca_input_path, "r") as orca_input:
            line = orca_input.readline()
            key = None
            while line and not line.startswith('  Spin-Determinant CI Printing'):
                line = orca_input.readline()
                if line.startswith('   Internal'):
                    self.internal = int(line.split()[5])
                if line.startswith('   Active'):
                    self.active = int(line.split()[5])
            if not line:
                raise ORCASectionNotFound('Spin-Determinant CI Printing')
            line = orca_input.readline()
            while line and not line.startswith('DENSITY MATRIX'):
                if line.startswith('CFG['):
                    key = line.split()[0][4:-1]
                    determinants[key] = []
                elif line.startswith('   ['):
                    det, val = line.split()
                    determinants[key].append({det[1:-1]: float(val)})
                line = orca_input.readline()

        with open(self.orca_input_path, "r") as orca_input:
            line = orca_input.readline()
            key = None
            while line and not line.startswith('  Extended CI Printing'):
                line = orca_input.readline()
            if not line:
                raise ORCASectionNotFound('Extended CI Printing')
            line = orca_input.readline()
            while line and not line.startswith('DENSITY MATRIX'):
                if line.startswith(' CFG['):
                    if determinants.get(line.split()[1]) == []:
                        key = line.split()[1]
                elif key and line.startswith(' \tCSF['):
                    determinants[key] = float(line.split()[1])
                    key = None
                line = orca_input.readline()

        self.spin_determinants = []
        # Sort spin-determinants by increasing occupation levels
        for k, v in sorted(determinants.items(), key=lambda x: ternary2decimal(x[0][::-1])):
            if isinstance(v, list):
                for v1 in v:
                    self.spin_determinants += [list(v1.items())[0]]
            else:
                self.spin_determinants += [(k, v)]

    def check_monotonic(self):
        """Check if energy of orbitals increases monotonically."""
        section_body = self.orca_section("ORBITAL ENERGIES")[3:]
        energy = None
        for line in section_body:
            split_line = line.split()
            if line.isspace():
                break
            elif len(split_line) == 4:
                if energy and float(split_line[2]) < energy:
                    print(line, energy)
                energy = float(split_line[2])

    @staticmethod
    def get_promotion_rules(spin_det_1, spin_det_2):
        """Creates promotions rules for two spin-determinants

        for '2u00' ->, '0u20'
        promotion rules is
            [(1, 1, 3), (2, 1, 3)]
        for rule (1, 1, 3)
            first is spin (1=up/2=down)
            second is number of orbital in active space electron promote from
            third is of orbital in active space electron promote to
        orbital numeration starting from 1.
        """
        # divide spin-determinants into spin parts
        spin_det_1_up = [e in ('2', 'u') for e in spin_det_1]
        spin_det_1_down = [e in ('2', 'd') for e in spin_det_1]

        spin_det_2_up = [e in ('2', 'u') for e in spin_det_2]
        spin_det_2_down = [e in ('2', 'd') for e in spin_det_2]

        # get difference between initial and final spin-det
        spin_det_up = [x - y for x, y in zip(spin_det_1_up, spin_det_2_up)]
        spin_det_down = [x - y for x, y in zip(spin_det_1_down, spin_det_2_down)]

        # get indexes of releasing and filling orbitals
        spin_det_1_up_indexes = [i+1 for i, x in enumerate(spin_det_up) if x == 1]
        spin_det_2_up_indexes = [i+1 for i, x in enumerate(spin_det_up) if x == -1]

        spin_det_1_down_indexes = [i+1 for i, x in enumerate(spin_det_down) if x == 1]
        spin_det_2_down_indexes = [i+1 for i, x in enumerate(spin_det_down) if x == -1]

        return (
            [(1, fr, to) for fr, to in zip(spin_det_1_up_indexes, spin_det_2_up_indexes)] +
            [(2, fr, to) for fr, to in zip(spin_det_1_down_indexes, spin_det_2_down_indexes)]
        )

    @property
    def ground_state(self):
        """Lowest occupation spin-determinant"""
        return self.spin_determinants[0][0]

    def correlation(self):
        """
        :returns: MDET section of correlation.data file
        """
        with open('correlation.data', 'w') as output_file:
            print('START MDET', file=output_file)
            print('Title', file=output_file)
            print(' multideterminant WFN %s\n' % self.title, file=output_file)

            print('MD', file=output_file)
            print('  %i' % len(self.spin_determinants), file=output_file)

            for i, (_, weight) in enumerate(self.spin_determinants):
                print(' %9.6f  %i %i' % (weight, i+1, int(i > 0)), file=output_file)

            for i, (spin_det, _) in enumerate(self.spin_determinants):
                if self.ground_state != spin_det:
                    for s, f, t in self.get_promotion_rules(self.ground_state, spin_det):
                        print('  DET %i %i PR %i 1 %i 1' % (i+1, s, f + self.internal, t + self.internal), file=output_file)
            print('END MDET', file=output_file)


def main():
    parser = argparse.ArgumentParser(
        description="This script converts multideterminant information to a CASINO correlation.data file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'code', type=int, help=(
            "number corresponding to the quantum chemistry code used to produce this MOLDEN file:\n"
            "3 -- ORCA 4.1\n"
            "7 -- QChem\n"
        )
    )
    parser.add_argument('input_file', type=str, help="path to output file")
    # truncation order for multideterminant extension
    parser.add_argument('--truncate', type=int, default=0, nargs='?', help="truncation order")
    parser.add_argument('--tolerance', type=float, default=0.01, nargs='?', help="min amplitude weight")
    args = parser.parse_args()

    if args.code == 0:
        Default().correlation()

    if args.code == 3:
        Orca(args.input_file).correlation()

    if args.code == 7:
        QChem(args.input_file, args.truncate, args.tolerance).correlation()


if __name__ == "__main__":
    main()
