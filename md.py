#!/usr/bin/env python2

import argparse
import copy
from yaml import dump, Dumper
from yaml.events import SequenceEndEvent
from collections import OrderedDict

parser = argparse.ArgumentParser(
    description="This script creates MULTIDETERMINANT INFORMATION section in gwfn.data file.",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument('prefix', type=str, default='mol', nargs='?', help="prifix of ORCA task")

args = parser.parse_args()


def get_determinants(prefix):
    determinants = OrderedDict()
    internal = 0
    active = 0

    with open(prefix + ".out", "r") as gwfn:
        line = gwfn.readline()
        key = None
        while line and not line.startswith('  Spin-Determinant CI Printing'):
            line = gwfn.readline()
            if line.startswith('   Internal'):
                internal = int(line.split()[5])
            if line.startswith('   Active'):
                active = int(line.split()[5])
        if not line:
            raise Exception
        line = gwfn.readline()
        while line and not line.startswith('DENSITY MATRIX'):
            if line.startswith('CFG['):
                key = line.split()[0][4:-1]
                determinants[key] = []
            elif line.startswith('   ['):
                det, val = line.split()
                determinants[key].append({det[1:-1]: float(val)})
            line = gwfn.readline()

    with open(prefix + ".out", "r") as gwfn:
        line = gwfn.readline()
        key = None
        while line and not line.startswith('  Extended CI Printing'):
            line = gwfn.readline()
        if not line:
            raise Exception
        line = gwfn.readline()
        while line and not line.startswith('DENSITY MATRIX'):
            if line.startswith(' CFG['):
                if determinants.get(line.split()[1]) == []:
                    key = line.split()[1]
            elif key and line.startswith(' \tCSF['):
                determinants[key] = float(line.split()[1])
                key = None
            line = gwfn.readline()

    result = []
    for k, v in determinants.items():
        if isinstance(v, list):
            for v1 in v:
                result += [v1.items()[0]]
        else:
            result += [(k, v)]

    return result, internal


def get_promote(c_1, c_2):

    c_1_up = [e in ('2', 'u') for e in c_1]
    c_1_down = [e in ('2', 'd') for e in c_1]

    c_2_up = [e in ('2', 'u') for e in c_2]
    c_2_down = [e in ('2', 'd') for e in c_2]

    c_up = [x - y for x, y in zip(c_1_up, c_2_up)]
    c_down = [x - y for x, y in zip(c_1_down, c_2_down)]

    c_1_up_indexes = [i+1 for i, x in enumerate(c_up) if x == 1]
    c_2_up_indexes = [i+1 for i, x in enumerate(c_up) if x == -1]

    c_1_down_indexes = [i+1 for i, x in enumerate(c_down) if x == 1]
    c_2_down_indexes = [i+1 for i, x in enumerate(c_down) if x == -1]

    return [(1, f, t) for f, t in zip(c_1_up_indexes, c_2_up_indexes)] + [(2, f, t) for f, t in zip(c_1_down_indexes, c_2_down_indexes)]

dets, internal = get_determinants(args.prefix)

first_c = dets[0][0]

print '  MD'
print '  %i' % len(dets)
for _, b in dets:
    print ' %9.6f' % b

for i, (a, _) in enumerate(dets):
    if first_c != a:
        for s, f, t in get_promote(first_c, a):
            print '  DET %i %i PR %i 1 %i 1' % (i+1, s, f + internal, t + internal)
