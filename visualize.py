#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
https://bitbucket.org/gutenkunstlab/dadi
"""
import os
import re
import json

import dadi

import matplotlib.pyplot as plt
import seaborn as sns

import run_dadi

sns.set_style('white')
#  darkgrid, whitegrid, dark, white, ticks
#sns.set_palette('YlGnBu')
#sns.despine()


#########1#########2#########3#########4#########5#########6#########7#########
# Visualization

def heatmap_sfs2d(fs, ax=None, vmax=None):
    if not vmax:
        vmax = fs.max()
    output = sns.heatmap(fs, mask=fs.mask,
                         vmin=1, vmax=vmax,
                         norm=sns.mplcol.LogNorm(vmin=1, vmax=vmax),
                         ax=ax, square=True, cmap='Spectral')  # 'YlOrRd')
    output.invert_yaxis()
    return output


def plot_dadi2d(fs_obs, model_output):
    sns.plt.clf()
    fs_exp = dadi.Inference.optimally_scaled_sfs(model_output, fs_obs)
    fs_res = dadi.Inference.Anscombe_Poisson_residual(fs_exp, fs_obs,
                                                      mask=True)
    vmax = max(fs_obs.max(), fs_exp.max())
    fig, grid = sns.plt.subplots(2, 2, figsize=(12, 12))
    ((ax_obs, ax_exp), (ax_res, ax_hist)) = grid
    #ax_cbar = fig.add_axes([0.05, 0.4, 0.01, 0.3])
    #ax_cbar_res = fig.add_axes([0.93, 0.4, 0.01, 0.3])
    heatmap_sfs2d(fs_obs, ax_obs, vmax)
    heatmap_sfs2d(fs_exp, ax_exp, vmax)
    sns.heatmap(fs_res, mask=fs_res.mask,
                vmin=-abs(fs_res).max(), vmax=abs(fs_res.max()),
                ax=ax_res, square=True, cmap='RdBu_r').invert_yaxis()
    sns.distplot(fs_res.compressed(), ax=ax_hist, color='gray',
                 kde=False, rug=True)
    ax_obs.set_title('Observation')
    ax_exp.set_title('Expectation')
    ax_res.set_title('Residuals')
    ax_hist.set_title('Residuals')
    ax_obs.set_xlabel('Chi')
    ax_exp.set_xlabel('Chi')
    ax_res.set_xlabel('Chi')
    ax_obs.set_ylabel('Flo')
    ax_exp.set_ylabel('Flo')
    ax_res.set_ylabel('Flo')
    return fig


def save_png_sfs(fs_file):
    freq_spectrum = dadi.Spectrum.from_file(fs_file)
    fs_png = re.sub(r'\..+$', '.png', fs_file)
    dadi.Plotting.plot_single_2d_sfs(freq_spectrum, vmin=0.5)
    plt.savefig(fs_png)
    plt.close()


def plot_2d_comp_multinom(fs_exp, fs_obs):
    # Plot a comparison of the resulting fs with the data.
    import pylab
    pylab.figure(1)
    dadi.Plotting.plot_2d_comp_multinom(fs_exp, fs_obs, vmin=1, resid_range=3)


def growth(t, nu2b=0.0001, nu2f=0.01, T=100.0):
    return nu2b * (nu2f / nu2b) ** (t / T)


def test_growth():
    sns.plt.plot(range(100), [growth(t) for t in range(100)])
    sns.plt.draw()


#########1#########2#########3#########4#########5#########6#########7#########
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-m', '--mode',
                        choices=['full', 'time', 'minimal'], default='minimal')
    parser.add_argument('fsfile')
    parser.add_argument('paramfile')
    args = parser.parse_args()

    if args.mode == 'full':
        model = run_dadi.exponential_full_model
    elif args.mode == 'time':
        model = run_dadi.exponential_time_model
    elif args.mode == 'minimal':
        model = run_dadi.exponential_minimal_model
    extrap_log = dadi.Numerics.make_extrap_log_func(model)

    fs_obs = dadi.Spectrum.from_file(args.fsfile)
    with open(args.paramfile, 'r') as fin:
        param = json.load(fin)
    print(param)
    (root, ext) = os.path.splitext(args.paramfile)
    mu = float(root.split('_')[1])
    run_dadi.set_global(mu)
    prefix = 'popt-{}-{}'.format(root, args.mode)
    outfile = '{}_{:.2e}.png'.format(prefix, args.mutation)

    model_output = extrap_log(param, fs_obs.sample_sizes, run_dadi.pts_l)
    fig = plot_dadi2d(fs_obs, model_output)
    sns.plt.draw()
    fig.savefig(outfile)
    fig.clf()
