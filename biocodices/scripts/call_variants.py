#!/usr/bin/env python

import sys
from termcolor import colored
from os.path import expanduser
from biocodices import Sequencing


def call_variants_for_all_samples(base_dir):
    sequencing = Sequencing(expanduser(base_dir))
    print(biocodices_logo())
    print('Welcome to biocodices! Anlyzing reads for:')
    print(colored(sequencing, 'green'))
    print()

    for sample in sequencing.samples():
        sample.call_variants()
        print()

    print(colored('\nDone! Bless your heart.\n', 'green'))


def biocodices_logo():
    return (r"""
 _     _                     _ _
| |__ (_) ___   ___ ___   __| (_) ___ ___  ___
| '_ \| |/ _ \ / __/ _ \ / _` | |/ __/ _ \/ __|
| |_) | | (_) | (_| (_) | (_| | | (_|  __/\__ \
|_.__/|_|\___/ \___\___/ \__,_|_|\___\___||___/
""")


if __name__ == '__main__':
    args = sys.argv
    if len(args) != 2:
        print('\nUsage:')
        print('{} <root dir of sequencing>'.format(__file__))
        sys.exit()

    base_dir = args[1]
    call_variants_for_all_samples(base_dir)
