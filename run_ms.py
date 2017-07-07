#!/usr/bin/env python
import math
import os
import sys
import subprocess
import json

import run_dadi


def make_cmd(nrep, nu2b, nu2f):
    sample_sizes = [16, 16]
    L = 10000
    theta = run_dadi.theta_Florida
    T_abs = run_dadi.T_split
    N0 = run_dadi.N0
    T_ms = T_abs / (4 * N0)
    alpha = math.log(nu2f / nu2b) / T_ms
    strn = 'N0: {}, N2b: {}, N2f: {}\n'.format(N0, nu2b * N0, nu2f * N0)
    sys.stderr.write(strn)

    cmd = ['ms', sum(sample_sizes), nrep]
    cmd.extend(['-t', theta * L])
    cmd.extend(['-I', 2] + sample_sizes)
    cmd.extend(['-n', 2, nu2f])
    cmd.extend(['-g', 2, alpha])
    cmd.extend(['-ej', T_ms, 2, 1])
    cmd.extend(['-eg', T_ms, 2, 0])
    return [str(x) for x in cmd]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-r', '--nrep', type=int, default=1)
    parser.add_argument('infile')
    args = parser.parse_args()

    (root, ext) = os.path.splitext(args.infile)
    mu = float(root.split('_')[1])
    run_dadi.set_global(mu)
    with open(args.infile, 'r') as fin:
        params = json.load(fin)
    cmd = make_cmd(args.nrep, *params)
    sys.stderr.write(' '.join(cmd) + '\n')
    if not args.dry_run:
        subprocess.call(cmd)
