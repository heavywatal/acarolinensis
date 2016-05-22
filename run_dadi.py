#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
https://bitbucket.org/gutenkunstlab/dadi
"""
import os
import re
import itertools
import json

import dadi

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')
#  darkgrid, whitegrid, dark, white, ticks
#sns.set_palette('YlGnBu')
#sns.despine()

np.set_printoptions(linewidth=160)


#########1#########2#########3#########4#########5#########6#########7#########
# Data

def marginal_stats(fs, dimension=0):
    d = dict()
    margin = fs.marginalize([dimension])
    d['pi'] = margin.pi()
    d['S'] = margin.S()
    d['D'] = margin.Tajima_D()
    d['theta_W'] = margin.Watterson_theta()
    d['theta_L'] = margin.theta_L()
    return d


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


def save_png_seaborn(outfile, fs_obs, model, param):
    sns.plt.clf()
    fig = plot_dadi2d(fs_obs, model(param, fs_obs.sample_sizes, pts_l))
    sns.plt.draw()
    fig.savefig(outfile)
    fig.clf()


def save_png_sfs(fs_file):
    freq_spectrum = dadi.Spectrum.from_file(fs_file)
    fs_png = re.sub(r'\..+$', '.png', fs_file)
    dadi.Plotting.plot_single_2d_sfs(freq_spectrum, vmin=0.5)
    plt.savefig(fs_png)
    plt.close()


def plot_2d_comp_multinom(fs_obs, func_ex):
    # Plot a comparison of the resulting fs with the data.
    import pylab
    pylab.figure(1)
    fs_exp = func_ex(p_opt, fs_obs.sample_sizes, pts_l)
    dadi.Plotting.plot_2d_comp_multinom(fs_exp, fs_obs, vmin=1, resid_range=3)


def growth(t, nu2b=0.0001, nu2f=0.01, T=100.0):
    return nu2b * (nu2f / nu2b) ** (t / T)


def test_growth():
    sns.plt.plot(range(100), [growth(t) for t in range(100)])
    sns.plt.draw()


#########1#########2#########3#########4#########5#########6#########7#########
# Model

def exponential_full_model(params, fs_size, pts):
    """
    nu2b: The bottleneck size for pop2
    nu2f: The final size for pop2
    m12: The scaled migration rate from pop2 to pop1
    m21: The scaled migration rate from pop1 to pop2
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    (nu2b, nu2f, m12, m21, T) = params
    (n1, n2) = fs_size
    # Define the grid we'll use
    grid = dadi.Numerics.default_grid(pts)
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(grid)
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(grid, phi)

    def nu2_func(t):
        return nu2b * (nu2f / nu2b) ** (t / T)
    phi = dadi.Integration.two_pops(phi, grid, T, nu1=1, nu2=nu2_func,
                                    m12=m12, m21=m21)
    # Finally, calculate the spectrum.
    return dadi.Spectrum.from_phi(phi, (n1, n2), (grid, grid)).fold()


def exponential_time_model(params, fs_size, pts):
    (nu2b, nu2f, T) = params
    params = (nu2b, nu2f, 0, 0, T)
    return exponential_full_model(params, fs_size, pts)


#########1#########2#########3#########4#########5#########6#########7#########
# Exhaustive loops

def log_likelihood(fs_obs, func_ex, params):
    fs_model = func_ex(params, fs_obs.sample_sizes, pts_l)
    return dadi.Inference.ll_multinom(fs_model, fs_obs)


def exhaustive_loops(fs_obs, func_ex, params_grid):
    results = {}
    it = itertools.product(*params_grid)
    for params in it:
        print(params)
        results[params] = log_likelihood(fs_obs, func_ex, params)
    return results


def save_results(results, outfile):
    with open(outfile, 'w') as fout:
        json.dump(results.items(), fout)


#########1#########2#########3#########4#########5#########6#########7#########
# Parameters

# lizards LCNS: 5e-6 (Janes 2011)
# lizards synonymous: 5e-7 (Munoz 2013)
# lizards mtDNA: 6.5e-9 = 0.65% per MY (Macey)
# lizards nDNA genes: 2.1e-10 = 0.021% per MY (Tollis 2014 Genetica)
# avian: 1.23--2.21 e-9 per site per year (Nam 2010 Genome Biol)
# human: 2.5e-8 (Nachman 2000 Genetics)
# mammalian: 2.2e-9 per base per year (Kumar 2002 PNAS)
# Drosophila: 3.5e-9 per base generation (Keightley 2009 Genome Res)
# nematodes: 2.1e-8 per site per generation (Denver 2004 Nature)
# Arabidopsis: 7.1e-9 per base per generation (Ossowski 2010 Science)
# E. coli: 2.2e-10 per base per generation (Lee 2012 PNAS)

pts_l = [40, 50, 60]
T_split = 50.0
theta_Florida = 0.0014


def calc_N0(u):
    # theta = 4Nu
    return theta_Florida / (4 * u)


def make_bounds_full(N0):
    # parameters (nu2b, nu2f, m12, m21, T)
    lower_bound = [2 / N0, 0.0001, 1e-2, 1e-2, 50 / (2 * N0)]
    upper_bound = [1000 / N0, 1.0, 10, 10, 500 / (2 * N0)]
    init = [100 / N0, 0.01, 1, 1, 100 / (2 * N0)]
    return (lower_bound, upper_bound, init)


def make_bounds_time(N0):
    # parameters (nu2b, nu2f, T)
    lower_bound = [2 / N0, 0.0001, 50 / (2 * N0)]
    upper_bound = [1000 / N0, 1.0, 500 / (2 * N0)]
    init = [100 / N0, 0.01, 100 / (2 * N0)]
    return (lower_bound, upper_bound, init)


def make_grid(lower_bound, upper_bound, breaks=6):
    z = zip(lower_bound, upper_bound)
    print(z)
    return [np.logspace(np.log10(l), np.log10(u), breaks) for (l, u) in z]


def name_params(params, u):
    return {'u': u,
            'N0': calc_N0(u),
            'nu2b': params[0],
            'nu2f': params[1],
            'm12': params[2],
            'm21': params[3],
            'T': params[4]}


def translate(params):
    N0 = params['N0']
    params['nu2b'] *= N0
    params['nu2f'] *= N0
    params['m12'] *= N0 * 2
    params['m21'] *= N0 * 2
    params['T'] *= N0 * 2
    return params


#########1#########2#########3#########4#########5#########6#########7#########
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-e', '--exhaustive', action='store_true')
    parser.add_argument('-o', '--optimize', action='store_true')
    parser.add_argument('-f', '--full', action='store_true')
    parser.add_argument('-u', '--mutation', type=float, default=1e-8)
    parser.add_argument('-b', '--breaks', type=int, default=6)
    parser.add_argument('-l', '--load')
    parser.add_argument('infile')
    args = parser.parse_args()

    (root, ext) = os.path.splitext(args.infile)
    fs_obs = dadi.Spectrum.from_file(args.infile)

    # Make the extrapolating version of our demographic model function.
    if args.full:
        model = exponential_full_model
        fixed = [None, None, None, None, None]
        make_bounds = make_bounds_full
    else:
        model = exponential_time_model
        fixed = [None, None, None]
        make_bounds = make_bounds_time
    extrap_log = dadi.Numerics.make_extrap_log_func(model)

    if args.load:
        (base, ext) = os.path.splitext(args.load)
        u = float(base.split('_')[1])
        N0 = calc_N0(u)
        (lower_bound, upper_bound, p0) = make_bounds(N0)
        with open(args.load, 'r') as fin:
            p0 = json.load(fin)
    else:
        N0 = calc_N0(args.mutation)
        (lower_bound, upper_bound, p0) = make_bounds(N0)

    T_rel = T_split / (2 * N0)
    print('theta={}, u={}, N0={}, T={}'.format(
        theta_Florida, args.mutation, N0, T_rel))
    print(lower_bound)
    print(upper_bound)

    if args.dry_run:
        exit()
    elif args.optimize:
        p_opt = dadi.Inference.optimize_log(p0, fs_obs, extrap_log, pts_l,
                                            lower_bound=lower_bound,
                                            upper_bound=upper_bound,
                                            fixed_params=fixed,
                                            epsilon=2e-3,
                                            verbose=1, maxiter=len(p0) * 200)
        print('log(lik): ' + str(log_likelihood(fs_obs, extrap_log, p_opt)))
        p_opt = p_opt.tolist()
        print(p_opt)
        # params = name_params(p_opt, args.mutation)
        # print(translate(params))
        if args.full:
            prefix = 'popt-' + root + '-full'
        else:
            prefix = 'popt-' + root + '-time'
        if args.load:
            prefix += '-loadpexh'
        with open(prefix + '.json', 'w') as fout:
            json.dump(p_opt, fout)
        outfile = prefix + '.png'
        save_png_seaborn(outfile, fs_obs, extrap_log, p_opt)
    elif args.exhaustive:
        params_grid = make_grid(lower_bound, upper_bound, args.breaks)
        print(params_grid)
        if args.full:
            prefix = root + '-full'
        else:
            prefix = root + '-time'
        outfile = '{}_{:.2e}.json'.format(prefix, args.mutation)
        print(outfile)
        results = exhaustive_loops(fs_obs, extrap_log, params_grid)
        save_results(results, outfile)
    else:
        print(marginal_stats(fs_obs, 0))
        print(marginal_stats(fs_obs, 1))
        if not args.load:
            exit()
        p_opt = p0
        fs_model = extrap_log(p_opt, fs_obs.sample_sizes, pts_l)
        print(marginal_stats(fs_model, 0))
        print(marginal_stats(fs_model, 1))
        fs_exp = dadi.Inference.optimally_scaled_sfs(fs_model, fs_obs)
        print(marginal_stats(fs_exp, 0))
        print(marginal_stats(fs_exp, 1))
        outfile = root + '_test.png'
        save_png_seaborn(outfile, fs_obs, extrap_log, p_opt)
