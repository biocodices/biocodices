#!/usr/bin/env python

import sys
from os.path import basename
from collections import defaultdict
from biocodices.helpers.language import plural


def count_variants_and_samples(vcf_path):
    with open(vcf_path, 'r') as vcf_file:
        lines =  [line.strip() for line in vcf_file.readlines()]

    variants = [line for line in lines if not line.startswith('#')]

    filters = defaultdict(int)
    for line in variants:
        filter_status = line.split('\t')[6]
        filters[filter_status] += 1

    comments = [line for line in lines if line.startswith('#')]

    header = [line for line in comments if line.startswith('#CHROM')][0]
    samples = header.split('\t')[9:]

    print('{}, filters: {} \n{}'.format(plural('variant', len(variants)),
                                        dict(filters),
                                        plural('sample', len(samples))))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('\nUsage:')
        print('{} <vcf file>\n'.format(basename(__file__)))
        sys.exit()

    vcf_path = sys.argv[1]
    count_variants_and_samples(vcf_path)
