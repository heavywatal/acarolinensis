#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import json

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import run_dadi

sns.set_style('white')
#########1#########2#########3#########4#########5#########6#########7#########


def make_df(infile):
    with open(infile, 'r') as fin:
        results = json.load(fin)
    tr = [run_dadi.translate(run_dadi.name_params(k)) for k, v in results]
    df = pd.DataFrame.from_records(tr)
    df['loglik'] = [v for k, v in results]
    return df


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

    table = make_df(args.infile)
    table = table[table.loglik > -5e5]
    print(table)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(10, 1e6)
    # cmap = sns.cubehelix_palette(start=1.5, rot=0.3, as_cmap=True)
    cmap = plt.cm.get_cmap('YlGnBu')
    kws = {'c': table['loglik'], 's': 100, 'cmap': cmap}
    reg = sns.regplot('nu2b', 'nu2f', data=table, ax=ax, fit_reg=False,
                      marker='s', scatter_kws=kws)
    # print(out.get_children())
    path_collection = reg.get_children()[0]
    plt.colorbar(mappable=path_collection)
    # plt.show()  # opens a window
    fig.savefig('loglik_{}.png'.format(root))
