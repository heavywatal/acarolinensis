#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
https://bitbucket.org/gutenkunstlab/dadi
"""
import os
import itertools
import json

import dadi

import numpy as np
np.set_printoptions(linewidth=160)


# #######1#########2#########3#########4#########5#########6#########7#########
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


# #######1#########2#########3#########4#########5#########6#########7#########
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
        return nu2b * ((nu2f / nu2b) ** (t / T))
    phi = dadi.Integration.two_pops(phi, grid, T, nu1=1, nu2=nu2_func,
                                    m12=m12, m21=m21)
    # Finally, calculate the spectrum.
    return dadi.Spectrum.from_phi(phi, (n1, n2), (grid, grid)).fold()


def exponential_time_model(params, fs_size, pts):
    (nu2b, nu2f, T) = params
    params = (nu2b, nu2f, 0, 0, T)
    return exponential_full_model(params, fs_size, pts)


def exponential_minimal_model(params, fs_size, pts):
    (nu2b, nu2f) = params
    params = (nu2b, nu2f, 0, 0, T_rel)
    return exponential_full_model(params, fs_size, pts)


# #######1#########2#########3#########4#########5#########6#########7#########
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


# #######1#########2#########3#########4#########5#########6#########7#########
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
N0 = None
T_split = 50.0
T_rel = None
theta_Florida = 0.0014


def set_global(u):
    # theta = 4Nu
    global N0
    global T_rel
    N0 = theta_Florida / (4 * u)
    T_rel = T_split / (2.0 * N0)


def make_bounds_full():
    # parameters (nu2b, nu2f, m12, m21, T)
    lower_bound = [2 / N0, 0.0001, 1e-2, 1e-2, 50 / (2 * N0)]
    upper_bound = [1000 / N0, 1.0, 10, 10, 500 / (2 * N0)]
    init = [100 / N0, 0.01, 1, 1, 100 / (2 * N0)]
    return (lower_bound, upper_bound, init)


def make_bounds_time():
    # parameters (nu2b, nu2f, T)
    lower_bound = [2 / N0, 0.001, 50 / (2 * N0)]
    upper_bound = [1000 / N0, 1.0, 100 / (2 * N0)]
    init = [100 / N0, 0.01, 100 / (2 * N0)]
    return (lower_bound, upper_bound, init)


def make_bounds_minimal():
    # parameters (nu2b, nu2f, T)
    lower_bound = [2 / N0, 0.001]
    upper_bound = [1000 / N0, 1.0]
    init = [100 / N0, 0.01]
    return (lower_bound, upper_bound, init)


def make_grid(lower_bound, upper_bound, breaks):
    z = zip(lower_bound, upper_bound)
    print(z)
    return [np.logspace(np.log10(l), np.log10(u), breaks) for (l, u) in z]


def name_params(params):
    names = ('nu2b', 'nu2f', 'm12', 'm21', 'T')
    if len(params) <= 3:
        names = names[:2]
    if len(params) == 3:
        names.append('T')
    return {k: v for k, v in zip(names, params)}


def translate(params):
    params['nu2b'] *= N0
    params['nu2f'] *= N0
    params['m12'] = params.get('m12', 0) * N0 * 2
    params['m21'] = params.get('m21', 0) * N0 * 2
    params['T'] = params.get('T', T_rel) * N0 * 2
    return params


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-e', '--exhaustive', action='store_true')
    parser.add_argument('-o', '--optimize', action='store_true')
    parser.add_argument('-m', '--mode',
                        choices=['full', 'time', 'minimal'], default='minimal')
    parser.add_argument('-u', '--mutation', type=float, default=1e-8)
    parser.add_argument('-b', '--breaks', type=int, default=20)
    parser.add_argument('-l', '--load')
    parser.add_argument('infile')
    args = parser.parse_args()

    (root, ext) = os.path.splitext(args.infile)
    fs_obs = dadi.Spectrum.from_file(args.infile)

    # Make the extrapolating version of our demographic model function.
    if args.mode == 'full':
        model = exponential_full_model
        fixed = [None, None, None, None, None]
        make_bounds = make_bounds_full
    elif args.mode == 'time':
        model = exponential_time_model
        fixed = [None, None, None]
        make_bounds = make_bounds_time
    elif args.mode == 'minimal':
        model = exponential_minimal_model
        fixed = [None, None]
        make_bounds = make_bounds_minimal
    extrap_log = dadi.Numerics.make_extrap_log_func(model)

    if args.load:
        (base, ext) = os.path.splitext(args.load)
        args.mutation = float(base.split('_')[1])
        set_global(args.mutation)
        (lower_bound, upper_bound, p0) = make_bounds()
        with open(args.load, 'r') as fin:
            p0 = json.load(fin)
    else:
        set_global(args.mutation)
        (lower_bound, upper_bound, p0) = make_bounds()

    print('theta={}, u={}, N0={}, T={}'.format(
        theta_Florida, args.mutation, N0, T_rel))
    print(lower_bound)
    print(upper_bound)

    if args.optimize:
        p_opt = dadi.Inference.optimize_log(p0, fs_obs, extrap_log, pts_l,
                                            lower_bound=lower_bound,
                                            upper_bound=upper_bound,
                                            fixed_params=fixed,
                                            epsilon=1e-3,
                                            verbose=1, maxiter=len(p0) * 200)
        print('log(lik): ' + str(log_likelihood(fs_obs, extrap_log, p_opt)))
        p_opt = p_opt.tolist()
        print(p_opt)
        # params = name_params(p_opt, args.mutation)
        # print(translate(params))
        prefix = 'popt-{}-{}'.format(root, args.mode)
        if args.load:
            prefix += '-loadpexh'
        outfile = '{}_{:.2e}.json'.format(prefix, args.mutation)
        print('Writing ' + outfile)
        with open(outfile, 'w') as fout:
            json.dump(p_opt, fout)
    elif args.exhaustive:
        params_grid = make_grid(lower_bound, upper_bound, args.breaks)
        print(params_grid)
        prefix = '{}-{}'.format(root, args.mode)
        outfile = '{}_{:.2e}.json'.format(prefix, args.mutation)
        print('Writing ' + outfile)
        if args.dry_run:
            exit()
        results = exhaustive_loops(fs_obs, extrap_log, params_grid)
        save_results(results, outfile)
    else:
        print('Note: add -o or -e to run dadi')
        print(fs_obs.sample_sizes)
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
