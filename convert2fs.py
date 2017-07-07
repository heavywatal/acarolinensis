#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import re
import dadi


def get_sample_size(data_dict):
    result = dict()
    for locus in data_dict.values():
        for pop, freqs in locus['calls'].items():
            result[pop] = max(sum(freqs), result.get(pop, 0))
    return result


def convert2fs(snp_file):
    """Parse the data file to generate the data dictionary
    """
    fs_file = re.sub(r'\..+$', '.fs', snp_file)
    data_dict = dadi.Misc.make_data_dict(snp_file)
    print('{} loci'.format(len(data_dict)))
    counts = get_sample_size(data_dict)
    pop_ids = sorted(counts.keys(), reverse=True)
    print(pop_ids)
    assert pop_ids == ['Flo', 'Chi']
    sample_sizes = [counts[x] for x in pop_ids]
    print(sample_sizes)
    freq_spectrum = dadi.Spectrum.from_data_dict(
        data_dict, pop_ids, sample_sizes, polarized=False)
    # polarized=False means folded=True
    print('writing: ' + fs_file)
    freq_spectrum.to_file(fs_file)
    return freq_spectrum


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    args = parser.parse_args()
    convert2fs(args.infile)
