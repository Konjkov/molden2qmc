#!/usr/bin/env python3
import re
import argparse
from functools import reduce
from collections import OrderedDict
from decimal import Decimal


class SectionNotFound(Exception):
    """Section not found in MOLDEN file."""

    def __init__(self, section_name):
        self.section_name = section_name

    def __str__(self):
        return repr(self.section_name)


class QChemSectionNotFound(SectionNotFound):
    """Section not found in QChem Output file."""


class ORCASectionNotFound(SectionNotFound):
    """Section not found in ORCA Output file."""


class Default:

    title = "Single determinant."

    def __init__(self, *args, **kwargs):
        """Single determinant"""
        self.determinants = [('', 1.0)]

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

    def ground_state(self, det):
        """Lowest occupation spin-determinant."""
        count_2 = det.count('2')
        count_u = det.count('u')
        count_d = det.count('d')
        count_0 = det.count('0')
        return (
            '2' * (count_2 + min(count_u, count_d)) +
            'u' * (count_u - count_d) +
            'd' * (count_d - count_u) +
            '0' * (count_0 + min(count_u, count_d))
        )

    def correlation(self):
        """
        :returns: empty MDET section of correlation.data file
        """
        with open('correlation.data', 'w') as output_file:
            print('START MDET', file=output_file)
            print('Title', file=output_file)
            print(' multideterminant WFN {}\n'.format(self.title), file=output_file)
            print('MD', file=output_file)
            print('  {}'.format(1), file=output_file)
            print(' {: .9f}  {} {}'.format(1, 1, 0), file=output_file)
            print('END MDET', file=output_file)


class PSI4(Default):
    """
    Psi4
    """
    title = "generated from Psi4 output data."

    def __init__(self, input_path):
        """Initialise multi-determinant support."""
        super().__init__(input_path)
        self.occupied = {'alpha': 0, 'beta': 0}
        self.active = {'alpha': 0, 'beta': 0}
        self.input_path = input_path
        self.parse_output()

    def parse_output(self):
        """Retrive from Psi4 output:

           The 20 most important determinants:

            *   1    0.967659  (    0,    0)  2AX 3AA
            *   2    0.101346  (    2,    2)  3AA 4AX
            *   3    0.101346  (    4,    3)  3AA 5AX
            *   4    0.076684  (    4,    5)  3AA 5AA 7AB
            *   5    0.076684  (    2,    6)  3AA 4AA 8AB
            *   6    0.076581  (   11,    3)  3AA 5AB 7AA
            *   7    0.076581  (   16,    2)  3AA 4AB 8AA
            *   8   -0.071317  (   21,    1)  2AA 3AB 9AA
            *   9    0.058707  (   11,    5)  3AA 7AX
            *  10    0.058707  (   16,    6)  3AA 8AX
            *  11    0.041218  (    6,    7)  2AA 6AA 9AB
            *  12   -0.038805  (   22,    0)  2AB 3AA 9AA
            *  13   -0.034382  (    7,    1)  3AX 6AA
            *  14   -0.032511  (    0,    7)  2AA 3AA 9AB
            *  15    0.030375  (   22,    7)  3AA 9AX
            *  16    0.024633  (   21,    4)  2AA 6AB 9AA
            *  17    0.017823  (    7,    4)  3AA 6AX
            *  18   -0.016585  (   25,    0)  2AB 6AA 9AA
            *  19   -0.005510  (   25,    7)  6AA 9AX
            *  20    0.002691  (    6,    0)  2AX 6AA

        :return: list of spin-determinants, number of internal orbitals
        """
        with open(self.input_path, "r") as psi4_input:
            line = psi4_input.readline()
            key = None
            while line and not line.startswith('   The '):
                line = psi4_input.readline()
                if line.startswith('         Frozen DOCC'):
                    docc = int(line.split()[3])
                if line.startswith('              Active'):
                    active = int(line.split()[2])
            if not line:
                return
            self.determinants = []
            while line:
                line = psi4_input.readline()
                if line.startswith('    *'):
                    line = line.replace('(', '').replace(')', '').replace(',', '')
                    weight = Decimal(line.split()[2])
                    determinant_tmpl = ['2'] * docc + ['0'] * active
                    for orb in line.split()[5:]:
                        if orb[-2:] == 'AX':
                            determinant_tmpl[int(orb[:-2])-1] = '2'
                        elif orb[-2:] == 'AA':
                            determinant_tmpl[int(orb[:-2])-1] = 'u'
                        elif orb[-2:] == 'AB':
                            determinant_tmpl[int(orb[:-2])-1] = 'd'
                    self.determinants.append((''.join(determinant_tmpl), weight))

    def correlation(self):
        """
        :returns: MDET section of correlation.data file
        """
        with open('correlation.data', 'w') as output_file:
            print('START MDET', file=output_file)
            print('Title', file=output_file)
            print(' multideterminant WFN %s\n' % self.title, file=output_file)

            print('MD', file=output_file)
            print('  %i' % (len(self.determinants)), file=output_file)
            opt_group_number = 0
            prev_weight = None
            for i, (_, weight) in enumerate(self.determinants):
                if prev_weight != weight:
                    opt_group_number += 1
                print(' %9.6f  %i %i' % (weight, opt_group_number, int(i > 0)), file=output_file)
                prev_weight = weight
            for i, (spin_det, _) in enumerate(self.determinants):
                if self.ground_state(spin_det) != spin_det:
                    for s, f, t in self.get_promotion_rules(self.ground_state(spin_det), spin_det):
                        print('  DET %i %i PR %i 1 %i 1' % (i+1, s, f, t), file=output_file)
            print('END MDET', file=output_file)


class Orca(Default):
    """
    ORCA 4.X
    """
    title = "generated from Orca output data"

    def __init__(self, input_path):
        """Initialise multi-determinant support."""
        super().__init__(input_path)
        # ORCA output precision
        self.tolerance = Decimal('0.000000001')
        self.internal = 0  # CASSCF internal orbitals
        self.active = 0  # CASSCF active orbitals
        self.input_path = input_path
        self.parse_output()

    def orca_section(self, section_name):
        """
        :returns: content of named section
        """
        with open(self.input_path, "r") as orca_input:
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

        :return: list of spin-determinants, number of internal orbitals
        """

        determinants = dict()
        cfg = weight = None

        with open(self.input_path, "r") as orca_input:
            line = orca_input.readline()
            while line and not (line.startswith('  Spin-Determinant CI Printing') or line.startswith('  Extended CI Printing')):
                line = orca_input.readline()
                if line.startswith('   Internal'):
                    self.internal = int(line.split()[5])
                if line.startswith('   Active'):
                    self.active = int(line.split()[5])
            if not line:
                return
            line = orca_input.readline()
            while line and not line.startswith('DENSITY MATRIX'):
                # Group determinants by CFD weights
                if line.startswith('CFG['):
                    cfg, weight = line.split()
                    cfg = cfg[4:-1]
                    determinants[cfg] = dict(weight=Decimal(weight), spin_det=[])
                elif weight and line.startswith('   ['):
                    det, val = line.split()
                    spin_det = [det[1:-1], Decimal(val)]
                    if abs(spin_det[1]) > self.tolerance:
                        spin_det = self.set_promotion_parity(spin_det)
                        spin_det = self.set_determinant_parity(spin_det)
                        determinants[cfg]['spin_det'].append(spin_det)
                line = orca_input.readline()

        with open(self.input_path, "r") as orca_input:
            line = orca_input.readline()
            key = val = None
            while line and not line.startswith('  Extended CI Printing'):
                line = orca_input.readline()
            if not line:
                raise ORCASectionNotFound('Extended CI Printing')
            line = orca_input.readline()
            while line and not line.startswith('DENSITY MATRIX'):
                if line.startswith(' CFG['):
                    _, key, _, _ = line.split()
                    val = determinants.get(key)
                elif key and val and not val['spin_det'] and line.startswith(' \tCSF['):
                    _, spin_det = line.split()
                    val['spin_det'].append((key, Decimal(spin_det)))
                line = orca_input.readline()
        self.determinants = OrderedDict(sorted(determinants.items(), key=lambda x: abs(x[1]['weight']), reverse=True))

    def parity_shell(self, values):
        """Determine parity of list of integers.
        from http://www.dalkescientific.com/writings/diary/archive/2016/08/14/fragment_chiral_molecules.html
        """
        N = len(values)
        num_swaps = 0
        for i in range(N - 1):
            for j in range(i + 1, N):
                if values[i] > values[j]:
                    values[i], values[j] = values[j], values[i]
                    num_swaps += 1
        return num_swaps % 2

    def set_determinant_parity(self, det):
        """parity of determinant in CSF in ORCA
        '2u000000' -  1
        '0u200000' - -1
        '0u020000' - -1
        'du00u000' - -1
        'ud00u000' -  1
        'uu00d000' - -1
        """
        csf = []
        count_u = 0
        count_d = 0
        for x in det[0]:
            if x == 'u':
                csf.append(count_u * 2)
                count_u += 1
            elif x == 'd':
                csf.append(count_d * 2 + 1)
                count_d += 1
            elif x == '2':
                csf.append(count_u * 2)
                count_u += 1
                csf.append(count_d * 2 + 1)
                count_d += 1
        parity = self.parity_shell(csf)
        return det[0], (-1)**parity*det[1]

    def set_promotion_parity(self, det):
        """parity of promotion rule in CASINO"""
        count_2 = det[0].count('2')
        count_u = det[0].count('u')
        count_d = det[0].count('d')
        count_0 = det[0].count('0')
        up_det = list(range(count_2 + count_u)) + (count_0 + count_d) * [None]
        down_det = list(range(count_2 + count_d)) + (count_0 + count_u) * [None]
        rules = self.get_promotion_rules(self.ground_state(det[0]), det[0])
        for s, f, t in rules:
            if s == 1:
                up_det[t-1], up_det[f-1] = up_det[f-1], None
            elif s == 2:
                down_det[t-1], down_det[f-1] = down_det[f-1], None
        up_det = [x for x in up_det if x is not None]
        down_det = [x for x in down_det if x is not None]
        parity = self.parity_shell(up_det) + self.parity_shell(down_det)
        return det[0], (-1)**parity*det[1]

    def correlation(self):
        """
        :returns: MDET section of correlation.data file
        """
        with open('correlation.data', 'w') as output_file:
            print('START MDET', file=output_file)
            print('Title', file=output_file)
            print(' multideterminant WFN {}\n'.format(self.title), file=output_file)
            print('MD', file=output_file)
            print('  {}'.format(sum((len(self.determinants[det]['spin_det']) for det in self.determinants))), file=output_file)
            opt_group_number = 0
            prev_weight = None
            for i, det in enumerate(self.determinants):
                weight = self.determinants[det]['weight']
                if prev_weight != weight:
                    opt_group_number += 1
                for spin_det in self.determinants[det]['spin_det']:
                    print(' {: .9f}  {} {}'.format(spin_det[1], opt_group_number, int(i > 0)), file=output_file)
                prev_weight = weight

            i = 0
            for det in self.determinants:
                for spin_det in self.determinants[det]['spin_det']:
                    i += 1
                    for s, f, t in self.get_promotion_rules(self.ground_state(spin_det[0]), spin_det[0]):
                        print('  DET {} {} PR {} 1 {} 1'.format(i, s, f + self.internal, t + self.internal), file=output_file)
            print('END MDET', file=output_file)


class QChem(Default):
    """
    QChem 4.4
    """
    title = "generated from QChem output data."

    def __init__(self, output_path, excitation=0, amplitude=0):
        """Initialise multi-determinant support."""
        super().__init__(output_path, excitation=0, amplitude=0)
        self.output_path = output_path
        self.occupied = {'alpha': 0, 'beta': 0}
        self.active = {'alpha': 0, 'beta': 0}
        self.determinants = []
        self.parse_output()
        self.truncate(excitation, amplitude)

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
        with open(self.output_path, "r") as qchem_output:
            line = qchem_output.readline()
            while line and not line.startswith('       Value      i             j           ->   a             b'):
                m = re.search(electron_regexp, line)
                if m:
                    self.occupied['alpha'] = int(m.group('alpha'))
                    self.occupied['beta'] = int(m.group('beta'))
                line = qchem_output.readline()
            if not line:
                # single determinant output
                return
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

    def truncate(self, excitation=None, amplitude=None):
        """Leave only determinants with active space orbital
        number not greater then order."""
        determinants = []
        limit = excitation + self.occupied['alpha']
        for det in self.determinants:
            if excitation and excitation and (det['promotions'][0]['to'] > limit or det['promotions'][1]['to'] > limit):
                continue
            if amplitude and abs(det['weight']) < amplitude:
                continue
            determinants.append(det)
        self.determinants = determinants

    def correlation(self):
        """
        :returns: MDET section of correlation.data file
        """
        with open('correlation.data', 'w') as output_file:
            print('START MDET', file=output_file)
            print('Title', file=output_file)
            print(' multideterminant WFN %s\n' % self.title, file=output_file)

            print('MD', file=output_file)
            print('  %i' % (len(self.determinants) + 1), file=output_file)
            # first determinant
            print(' %9.6f    1    0' % 1.0, file=output_file)
            opt_group_number = 1
            prev_weight = None
            for i, det in enumerate(self.determinants):
                if prev_weight != det['weight']:
                    opt_group_number += 1
                print(' %9.6f    %i    1' % (det['weight'], opt_group_number), file=output_file)
                prev_weight = det['weight']

            for i, det in enumerate(self.determinants):
                for p in det['promotions']:
                    # starting from 2-nd determinant
                    print('  DET %i %i PR %i 1 %i 1' % (i + 2, p['spin'], p['from'], p['to']), file=output_file)
            print('END MDET', file=output_file)


def main():
    parser = argparse.ArgumentParser(
        description="This script converts multideterminant information to a CASINO correlation.data file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'code', type=int, help=(
            "number corresponding to the quantum chemistry code used to produce this MOLDEN file:\n"
            "1 -- PSI4\n"
            "3 -- ORCA 4.1\n"
            "7 -- QChem\n"
        )
    )
    parser.add_argument('input_file', type=str, help="path to output file")
    # truncation for multideterminant extension
    parser.add_argument('--excitation', type=int, default=0, nargs='?', help="max excitaion orbital number")
    parser.add_argument('--amplitude', type=float, default=0, nargs='?', help="min amplitude weight")
    args = parser.parse_args()

    if args.code == 0:
        Default(args.input_file).correlation()

    if args.code == 1:
        PSI4(args.input_file).correlation()

    if args.code == 3:
        Orca(args.input_file).correlation()

    if args.code == 7:
        QChem(args.input_file, args.excitation, args.amplitude).correlation()


if __name__ == "__main__":
    main()
