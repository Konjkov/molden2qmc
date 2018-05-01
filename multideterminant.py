#!/usr/bin/env python3
import re
import argparse


class SectionNotFound(Exception):
    """Section not found in MOLDEN file."""

    def __init__(self, section_name):
        self.section_name = section_name

    def __str__(self):
        return repr(self.section_name)


class QChemSectionNotFound(Exception):
    """Section not found in QChem Output file."""

    def __init__(self, section_name):
        self.section_name = section_name

    def __str__(self):
        return repr(self.section_name)


class QChem:
    """
    QChem 4.4
    """
    title = "generated from QChem output data."

    def __init__(self, output_file, order=0):
        """Initialise multi-determinant support."""

        self.output_file = output_file
        self.occupied = {'alpha': 0, 'beta': 0}
        self.active = {'alpha': 0, 'beta': 0}
        self.determinants = []
        self.parse_output()
        # self.truncate(order)

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

    def truncate(self, order, tol=0.01):
        """Leave only determinants with active space orbital
        number not greater then order."""
        determinants = []
        limit = order + self.occupied['alpha']
        for det in self.determinants:
            if det['promotions'][0]['to'] > limit or det['promotions'][1]['to'] > limit:
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


def main():
    parser = argparse.ArgumentParser(
        description="This script converts multideterminant information to a CASINO correlation.data file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'code', type=int, help=(
            "number corresponding to the quantum chemistry code used to produce this MOLDEN file:\n"
            "7 -- QChem\n"
        )
    )

    parser.add_argument('input_file', type=str, help="path to QChem output file")
    args = parser.parse_args()

    if args.code == 7:
        QChem(args.input_file).correlation()


if __name__ == "__main__":
    main()
