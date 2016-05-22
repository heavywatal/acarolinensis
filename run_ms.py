#!/usr/bin/env python
import math
import sys
import subprocess
import json


def old():
    nrep = 1000000
    N0 = 200000
    nu2b = 0.00014057  # *N0 = 28.11335961
    nu2f = 0.00114639  # *N0 = 229.27796913
    T = 0.00012709
    # T_abs = T * 2 * N0  #= 50.83583664
    # nu2b = nu2f * exp(-alpha * T)
    call(nrep, nu2b, nu2f, N0, T)


def call(nrep, nu2b, nu2f, N0, T, **kwargs):
    nsam = 32
    L = 10000
    theta = 0.008
    T_abs = T * 2 * N0  # = 50.83583664
    T_ms = T_abs / (4 * N0)
    alpha = math.log(nu2f / nu2b) / T_ms

    cmd = ['ms', nsam, nrep]
    cmd.extend(['-t', theta * L])
    cmd.extend(['-I', 2, 16, 16])
    cmd.extend(['-n', 2, nu2f])
    cmd.extend(['-g', 2, alpha])
    cmd.extend(['-ej', T_ms, 2, 1])
    cmd.extend(['-eg', T_ms, 2, 0])
    cmd = [str(x) for x in cmd]
    sys.stderr.write(' '.join(cmd) + '\n')
    subprocess.call(cmd)


#########1#########2#########3#########4#########5#########6#########7#########
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nrep', type=int, default=1)
    parser.add_argument('infile')
    args = parser.parse_args()

    with open(args.infile, 'r') as fin:
        params = json.load(fin)
    print(params)
    call(args.nrep, **params)
