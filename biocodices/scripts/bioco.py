#!/usr/bin/env python

import argparse
import yaml
import sys
from termcolor import colored
from os.path import expanduser, join, dirname


import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# This prevents matplotlib raising an exception when running biocodices on a =
# remote server with no X. This line has to be executed before importing pyplot

from biocodices import Cohort, software_name
from biocodices.helpers.general import touch_all_the_logs, biocodices_logo
from biocodices.helpers import Stopwatch


def call_variants_for_all_samples(args):
    print(biocodices_logo())
    print('Welcome to {}! Anlyzing reads for...'.format(software_name))
    cohort = Cohort(expanduser(args.seq_dir))

    if args.samples:
        sample_ids = args.samples.split(',')
        for sample_id in sample_ids:
            if sample_id not in [sample.id for sample in cohort.samples]:
                msg = '{} not found in this cohort.'
                print(msg.format(sample_id))
                sys.exit()

        cohort.samples = [sample for sample in cohort.samples
                          if sample.id in sample_ids]

    if args.skip_samples:
        sample_ids = args.skip_samples.split(',')
        for sample_id in sample_ids:
            if sample_id not in [sample.id for sample in cohort.samples]:
                msg = '{} not found in this cohort.'
                print(msg.format(sample_id))
                sys.exit()

        cohort.samples = [sample for sample in cohort.samples
                          if sample.id not in sample_ids]

    print(colored(cohort, 'green'))
    print('\nYou can follow the details of the process with:')
    # other option: `tail -n0 -f {}/{{*/,}}*.log`
    print('`tail -n0 -f {}/{{*/,}}*.log`\n'.format(cohort.results_dir))

    dir_list = [sample.results_dir for sample in cohort.samples]
    touch_all_the_logs(cohort.dir, dir_list)


    stopwatch = Stopwatch().start()

    cohort.call_variants(trim_reads=args.trim_reads,
                         align_reads=args.align_reads,
                         create_vcfs=args.create_vcfs,
                         joint_genotyping=args.joint_genotyping,
                         hard_filtering=args.hard_filtering)

    stopwatch.stop()
    runtime = stopwatch.nice_total_run_time
    print('\nThe whole process took {}.'.format(runtime))
    print('Done! Bless your heart.\n')


if __name__ == '__main__':
    cli_yaml_path = join(dirname(__file__), 'cli_arguments.yml')
    with open(cli_yaml_path, 'r') as cli_yaml_file:
        cli = yaml.load(cli_yaml_file)

    parser = argparse.ArgumentParser(description=cli['description'])
    for flag, options in cli['args'].items():
        long_name, action, required, default = options
        parser.add_argument('-' + flag, '--' + long_name, action=action,
                            required=required, help=cli['help_texts'][flag],
                            default=default)

    args = parser.parse_args()
    if args.complete_pipeline:
        for attr in ['trim_reads', 'align_reads', 'create_vcfs',
                     'joint_genotyping', 'hard_filtering']:
            setattr(args, attr, True)

    call_variants_for_all_samples(args)
