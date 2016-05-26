#!/usr/bin/env python

import argparse
import yaml
from termcolor import colored
from os.path import expanduser, join, dirname

from biocodices import Cohort, software_name
from biocodices.helpers.general import touch_all_the_logs, biocodices_logo


def call_variants_for_all_samples(args):
    cohort = Cohort(expanduser(args.seq_dir))

    print(biocodices_logo())
    print('Welcome to {}! Anlyzing reads for:'.format(software_name))
    print(colored(cohort, 'green'))
    print('\nYou can follow the details of the process with:')
    # other option: `tail -n0 -f {}/{{*/,}}*.log`
    print('`tail -n0 -f {}/*/*.log`\n'.format(cohort.results_dir))

    dir_list = [sample.results_dir for sample in cohort.samples]
    touch_all_the_logs(cohort.dir, dir_list)

    cohort.call_variants(trim_reads=args.trim_reads,
                         align_reads=args.align_reads,
                         create_vcfs=args.create_vcfs,
                         joint_genotyping=args.joint_genotyping,
                         hard_filtering=args.hard_filtering)

    print('\nDone! Bless your heart.\n')


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
    call_variants_for_all_samples(args)
