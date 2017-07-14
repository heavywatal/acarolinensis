#!/usr/bin/env python
"""
.nested = list.files(pattern='stats_.*.tsv.gz') %>%
    {tibble(file=.)} %>%
    dplyr::mutate(data= purrr::map(file, read_tsv)) %>%
    print()

.flat = .nested %>%
   dplyr::mutate(file= str_replace_all(file, '^stats_|\\.tsv\\.gz$', '')) %>%
   tidyr::separate(file, c('b', 'f', 'j'), '_') %>%
   dplyr::mutate_at(vars(b, f, j), function(x) {
       str_replace(x, '^\\w', '') %>% readr::parse_double()
   }) %>%
   tidyr::unnest() %>%
   print()

.p = .flat %>% ggplot(aes(D_2, group=j, fill=j))+
   geom_histogram(bins=30, position='identity', alpha=0.4)+
   scale_fill_distiller(palette='Spectral', trans='log')+
   facet_grid(f ~ b)+
   theme_bw()
.p
ggsave('bottleneck.pdf', .p, width=6, height=6)
"""
import math
import gzip
import subprocess


def make_cmd(nrep, sample_sizes, b, f, j=50):
    L = 10000
    theta = 0.0014
    N0 = 35000
    T_ms = j / (4.0 * N0)
    alpha = math.log(f / b) / T_ms
    # strn = 'N0: {}, N2b: {}, N2f: {}\n'.format(N0, nu2b * N0, nu2f * N0)
    # sys.stderr.write(strn)

    cmd = ['ms', sum(sample_sizes), nrep]
    cmd.extend(['-t', theta * L])
    cmd.extend(['-I', 2] + sample_sizes)
    cmd.extend(['-n', 2, f])
    cmd.extend(['-g', 2, alpha])
    cmd.extend(['-ej', T_ms, 2, 1])
    cmd.extend(['-eg', T_ms, 2, 0])
    return [str(x) for x in cmd]


def make_params():
    for j in [50, 500, 5000]:
        for b in [5e-5, 5e-4, 5e-3]:
            for f in [0.005, 0.05, 0.5]:
                yield {'b': b, 'f': f, 'j': j}


def make_label(params):
    v = []
    for key, val in sorted(params.items()):
        v.append('{}{}'.format(key, val))
    return '_'.join(v)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-r', '--nrep', type=int, default=1)
    parser.add_argument('-I', '--sample-sizes',
                        nargs=2, type=int, default=[16, 16])
    args = parser.parse_args()

    for params in make_params():
        outfile = 'stats_' + make_label(params) + '.tsv.gz'
        print(outfile)
        cmd = make_cmd(args.nrep, args.sample_sizes, **params)
        # sys.stderr.write(' '.join(cmd) + '\n')
        if not args.dry_run:
            ms = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out = subprocess.check_output(['sample_stats++'], stdin=ms.stdout)
            with gzip.open(outfile, 'w') as fout:
                fout.write(out)
