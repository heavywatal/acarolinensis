#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Read json | sort by likelihood | print head(n) and write MLE to a new json
"""
import os
import json

import run_dadi


def load_json(infile, n=20):
    with open(infile, 'r') as fin:
        results = json.load(fin)
    results = [[k, v] for k, v in results]
    return sorted(results, key=lambda x: x[1], reverse=True)[0:n]


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
        print([run_dadi.translate(run_dadi.name_params(key)), value])
    param, ll = results[0]
    outfile = 'pexh-' + root + '.json'
    print('Writing ' + outfile)
    with open(outfile, 'w') as fout:
        json.dump(param, fout)
