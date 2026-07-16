#!/usr/bin/env python3
"""Generate the GEMINAL section of parameters.casl for CASINO.

CASL is not standard YAML (bare constraint lines like 2^g_5,5=2^g_4,4
are invalid YAML scalars), so the section is generated as plain text.

A "shell" is a list of symmetry-equivalent orbital indices (3 for p,
5 for d); the LAST index of each shell is the reference component that
appears in Parameters, the rest are tied to it via Constraints.

Usage: make_geminal_casl.py [rundir] --channels p:2 d:2 --mirror
(rundir defaults to the current directory; without --channels an
HF-like section is generated: Geminal 1 with occupied MOs only)

Shells are derived from <rundir>/gwfn.data: MOs are classified by
(l, m-slot) from the coefficient blocks, grouped into levels of 2l+1,
components ordered by m-slot with slot 0 last (the reference).
Levels that are not rotationally closed (e.g. D2h-mixed a_g pairs)
are demoted to diagonal-only groups with a warning, since
component-wise off-diagonal ties would break the 1S invariant.
"""
import os
import sys

import numpy as np

__version__ = '1.1'

NAME_L = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

# shell types (s/sp/p/d/f... 1/2/3/4/5...) -> l
SHELL_MAP = {1: 0, 2: 1, 3: 1, 4: 2, 5: 3, 6: 4}

# CASINO premultiplies the d solid-harmonic constant factors into the
# orbital coefficients (but not for f and g); remove them before
# comparing component weights
D_PREMULTIPLIED_FACTOR = np.array((0.5, 3.0, 3.0, 3.0, 6.0))


def value(v, opt):
    return f'[ {v}, {"optimizable" if opt else "fixed"} ]'


def geminal_section(occupied, anchors, shells, seeds=None, mirror=False, diag_shells=None, diag_seeds=None):
    seeds = seeds or {}
    diag_shells = diag_shells or []
    diag_seeds = diag_seeds or {}
    ref = [s[-1] for s in shells]
    gems = (2, 3) if mirror else (2,)
    lines = [
        'GEMINAL:',
        '  Default g optimizability: fixed',
        '  Default c optimizability: fixed',
        '  Geminal 1:',
        '    Parameters:',
        f'      c: {value(1.0, False)}',
    ]
    for m in occupied:
        lines.append(f'      g_{m},{m}: {value(1.0, False)}')
    if not shells and not diag_shells:
        return '\n'.join(lines)
    lines += [
        '  Geminal 2:',
        '    Parameters:',
        f'      c: {value(1.0, False)}',
    ]
    for m in anchors:
        lines.append(f'      g_{m},{m}: {value(1.0, False)}')
    for i in range(len(shells)):
        lines.append(f'      g_{ref[i]},{ref[i]}: {value(seeds.get((i, i), 0.0), True)}')
    for i, s in enumerate(diag_shells):
        lines.append(f'      g_{s[-1]},{s[-1]}: {value(diag_seeds.get(i, 0.0), True)}')
    for i in range(len(shells)):
        for j in range(i + 1, len(shells)):
            if len(shells[i]) == len(shells[j]):
                lines.append(f'      g_{ref[i]},{ref[j]}: {value(seeds.get((i, j), 0.0), True)}')
    if mirror:
        lines += [
            '  Geminal 3:',
            '    Parameters:',
            f'      c: {value(-1.0, False)}',
        ]
    lines.append('  Constraints:')
    for s in shells + diag_shells:
        members = []
        for g in gems:
            members += [f'{g}^g_{m},{m}' for m in [s[-1]] + s[-2::-1]]
        lines.append('    ' + '='.join(members))
    for i in range(len(shells)):
        for j in range(i + 1, len(shells)):
            if len(shells[i]) != len(shells[j]):
                continue
            a, b = shells[i], shells[j]
            pairs = [(a[-1], b[-1])] + list(zip(a[-2::-1], b[-2::-1]))
            members = []
            for g in gems:
                for m, k in pairs:
                    members += [f'{g}^g_{m},{k}', f'{g}^g_{k},{m}']
            lines.append('    ' + '='.join(members))
    return '\n'.join(lines)


def read_gwfn(base_path):
    with open(os.path.join(base_path, 'gwfn.data')) as f:
        lines = f.read().splitlines()

    def block(header, n, cast):
        i = next(j for j, l in enumerate(lines) if l.startswith(header)) + 1
        data = []
        while len(data) < n:
            line = lines[i]
            if not line.startswith('--'):
                if cast is float:
                    # floats are written in 20-char fields with no separator
                    data += [line[k * 20 : (k + 1) * 20] for k in range(len(line) // 20)]
                else:
                    data += line.split()
            i += 1
        return [cast(x) for x in data[:n]]

    nelec = block('Number of electrons', 1, int)[0]
    nshell = block('Number of shells', 1, int)[0]
    nbf = block('Number of basis functions', 1, int)[0]
    shell_moments = [SHELL_MAP[t] for t in block('Code for shell types', nshell, int)]
    mo = block('ORBITAL COEFFICIENTS', nbf * nbf, float)
    mo_up = np.array(mo).reshape((nbf, nbf))
    p = 0
    for l in shell_moments:
        if l == 2:
            mo_up[:, p : p + 2 * l + 1] /= D_PREMULTIPLIED_FACTOR
        p += 2 * l + 1
    return nelec, shell_moments, mo_up


def mo_levels(base_path):
    nelec, shell_moments, mo_up = read_gwfn(base_path)
    lm = np.array([(l, m) for l in shell_moments for m in range(2 * l + 1)])
    neu = nelec // 2
    by_l = {}
    for i, c in enumerate(mo_up):
        c2 = c ** 2
        l = max(set(lm[:, 0]), key=lambda t: c2[lm[:, 0] == t].sum())
        w = np.array([c2[(lm[:, 0] == l) & (lm[:, 1] == s)].sum() for s in range(2 * l + 1)])
        slot = int(np.argmax(w))
        pure = w[slot] / w.sum() > 0.98
        by_l.setdefault(int(l), []).append((i + 1, slot, pure))
    levels = {}
    for l, mos in by_l.items():
        n = 2 * l + 1
        levels[l] = []
        for k in range(0, len(mos), n):
            chunk = mos[k:k + n]
            closed = all(p for _, _, p in chunk) and sorted(s for _, s, _ in chunk) == list(range(n))
            if closed:
                chunk = sorted(chunk, key=lambda t: t[1])
                chunk = chunk[1:] + chunk[:1]
            levels[l].append(([i for i, _, _ in chunk], closed))
    return neu, levels


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('rundir', nargs='?', default='.')
    parser.add_argument('--channels', nargs='+', default=[], metavar='l:n')
    parser.add_argument('--mirror', action='store_true')
    parser.add_argument('--anchors', nargs='+', type=int, default=None,
        help='override anchor MOs (e.g. Be: 1, since the 2s pair is replaced by the p block)')
    parser.add_argument('--seed', type=float, default=-0.05,
        help='seed for the first full-block diagonal (Be near-degeneracy: -0.19)')
    args = parser.parse_args()
    neu, levels = mo_levels(args.rundir)
    shells, diag_shells = [], []
    for channel in args.channels:
        name, n = channel.split(':')
        available = levels.get(NAME_L[name], [])
        if int(n) > len(available):
            sys.exit(f'error: requested {name}:{n} but gwfn.data has only {len(available)} {name} level(s)')
        for members, closed in available[:int(n)]:
            if closed:
                shells.append(members)
            else:
                diag_shells.append(members)
                print(f'warning: {name} level {members} is not rotationally closed, '
                      f'using diagonal-only ties', file=sys.stderr)
    seeds = {(i, i): v for i, v in enumerate([args.seed, -0.02][:len(shells)])}
    diag_seeds = {0: -0.02} if diag_shells else {}
    occupied = list(range(1, neu + 1))
    if args.anchors is None:
        correlated = {m for s in shells + diag_shells for m in s}
        anchors = [m for m in occupied if m not in correlated]
    else:
        anchors = args.anchors
    print('# generated by: ' + ' '.join([sys.argv[0].rsplit('/', 1)[-1]] + sys.argv[1:]) + f' (v{__version__})')
    print(geminal_section(occupied, anchors, shells, seeds, args.mirror, diag_shells, diag_seeds))
