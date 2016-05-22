#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
https://bitbucket.org/RyanGutenkunst/dadi
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


def save_png_seaborn(outfile, fs_obs, func_ex, param):
    sns.plt.clf()
    fig = plot_dadi2d(fs_obs, func_ex(p_opt, fs_obs.sample_sizes, pts_l))
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

def exponential_model(params, fs_size, pts):
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


#########1#########2#########3#########4#########5#########6#########7#########
# Exhaustive loops

def log_likelihood(fs_obs, func_ex, params):
    fs_model = func_ex(params, fs_obs.sample_sizes, pts_l)
    return dadi.Inference.ll_multinom(fs_model, fs_obs)


def exhaustive_loops(fs_obs, func_ex,
                     range_nu2b, range_nu2f, range_m12, range_m21, range_T):
    results = {}
    it = itertools.product(
        range_nu2b, range_nu2f, range_m12, range_m21, range_T)
    for params in it:
        print(params)
        results[params] = log_likelihood(fs_obs, func_ex, params)
    return results


def print_head(results, n=20):
    lst = sorted(results.items(), key=lambda x: x[1], reverse=True)
    for (key, value) in lst[1:n]:
        print([key, value])


def save_results(results, outfile):
    with open(outfile, 'w') as fout:
        json.dump(results.items(), fout)


def load_json(infile, n=20):
    with open(infile, 'r') as fin:
        results = json.load(fin)
    return sorted(results, key=lambda x: x[1], reverse=True)[1:n]


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


def make_bounds_pre(N0):
    # parameters (nu2b, nu2f, m12, m21, T)
    lower_bound = [2 / N0, 0.0001, 1e-6, 1e-6, 50 / (2 * N0)]
    upper_bound = [1000 / N0, 1.0, 1e-2, 1e-2, 500 / (2 * N0)]
    init = [100 / N0, 0.01, 1e-4, 1e-4, 100 / (2 * N0)]
    return (lower_bound, upper_bound, init)


def make_bounds(N0):
    # parameters (nu2b, nu2f, m12, m21, T)
    lower_bound = [2 / N0, 0.0001, 0, 0, 50 / (2 * N0)]
    upper_bound = [0.001, 1.0, 0, 0, 500 / (2 * N0)]
    return (lower_bound, upper_bound)

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
    parser.add_argument('-e', '--exhaustive', action='store_true')
    parser.add_argument('-o', '--optimize', action='store_true')
    parser.add_argument('-l', '--load', action='store_true')
    parser.add_argument('-u', '--mutation', type=float, default=1e-8)
    parser.add_argument('infile')
    args = parser.parse_args()

    (root, ext) = os.path.splitext(args.infile)

    N0 = calc_N0(args.mutation)
    T_rel = T_split / (2 * N0)
    (lower_bound, upper_bound, p0) = make_bounds_pre(N0)

    print('theta={}, u={}, N0={}, T={}'.format(
        theta_Florida, args.mutation, N0, T_rel))
    print(lower_bound)
    print(upper_bound)

    # Make the extrapolating version of our demographic model function.
    extrap_log = dadi.Numerics.make_extrap_log_func(exponential_model)

    if args.optimize:
        fs_obs = dadi.Spectrum.from_file(args.infile)
        fixed = [None, None, 0, 0, None]
        p_opt = dadi.Inference.optimize_log(p0, fs_obs, extrap_log, pts_l,
                                            lower_bound=lower_bound,
                                            upper_bound=upper_bound,
                                            fixed_params=fixed,
                                            epsilon=2e-3,
                                            verbose=1, maxiter=len(p0) * 200)
        print('log(lik): ' + str(log_likelihood(fs_obs, extrap_log, p_opt)))
        print(p_opt.tolist())
        print(translate(name_params(p_opt, args.mutation)))
    elif args.exhaustive:
        fs_obs = dadi.Spectrum.from_file(args.infile)
        params_grid = make_grid(lower_bound, upper_bound, 6)
        print(params_grid)
        outfile = '{}_{:.2e}.json'.format(root, args.mutation)
        print(outfile)
        results = exhaustive_loops(fs_obs, extrap_log, *params_grid)
        print_head(results)
        save_results(results, outfile)
    elif args.load:
        u = float(root.split('_')[1])
        N0 = calc_N0(u)
        (lower_bound, upper_bound) = make_bounds(N0)
        fs_obs = dadi.Spectrum.from_file(root.split('_')[0] + '.fs')
        results = load_json(args.infile, 10)
        for (key, value) in results:
            print([key, value])
        p0, ll = results[0]
        p_opt = dadi.Inference.optimize_log(p0, fs_obs, extrap_log, pts_l,
                                            lower_bound=lower_bound,
                                            upper_bound=upper_bound,
                                            fixed_params=fixed,
                                            epsilon=4e-3,
                                            verbose=1, maxiter=len(p0) * 200)
        params = name_params(p_opt, u)
        with open('popt-' + root + '.json', 'w') as fout:
            json.dump(params, fout)
        outfile = root + '.png'
        save_png_seaborn(outfile, fs_obs, extrap_log, p_opt)

    else:
        fs_obs = dadi.Spectrum.from_file(args.infile)
        #save_png_sfs(args.infile)
        print(marginal_stats(fs_obs, 0))
        print(marginal_stats(fs_obs, 1))
        p_opt = [0.00047331468067680183, 0.02294559868542178, 0.0, 0.0, 0.0007142857142857143]
        fs_exp = extrap_log(p_opt, fs_obs.sample_sizes, pts_l)
        print(marginal_stats(fs_exp, 0))
        print(marginal_stats(fs_exp, 1))
        outfile = root + '_test.png'
        save_png_seaborn(outfile, fs_obs, extrap_log, p_opt)

# log(lik): -124374.93640508532
# [0.00047331468067680183, 0.02294559868542178, 0.0, 0.0, 0.0007142857142857143]
# {'nu2b': 16.566013823688063, 'm21': 0.0, 'nu2f': 803.09595398976228, 'm12': 0.0, 'u': 1e-08, 'T': 50.0, 'N0': 35000.0}

# log(lik): -132401.558699
# [0.0009977610940156722, 0.13027017292711654, 0.0, 0.0, 0.002260726523303243]
# {'nu2b': 34.921638290548529, 'm21': 0.0, 'nu2f': 4559.4560524490789, 'm12': 0.0, 'u': 1e-08, 'T': 158.25085663122701, 'N0': 35000.0}
