#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import json

import run_dadi
#########1#########2#########3#########4#########5#########6#########7#########


def load_json(infile, n=20):
    with open(infile, 'r') as fin:
        results = json.load(fin)
    results = [[k, v] for k, v in results if k[0] < k[1]]
    return sorted(results, key=lambda x: x[1], reverse=True)[1:n]


#########1#########2#########3#########4#########5#########6#########7#########
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=30)
    parser.add_argument('infile')
    args = parser.parse_args()

    (root, ext) = os.path.splitext(args.infile)
    mu = float(root.split('_')[1])
    run_dadi.set_global(mu)

    results = load_json(args.infile, args.n)
    for (key, value) in results:
        print([run_dadi.translate(run_dadi.name_params(key, mu)), value])
        # print([key, value])
    param, ll = results[0]
    with open('pexh-' + root + '.json', 'w') as fout:
        json.dump(param, fout)
