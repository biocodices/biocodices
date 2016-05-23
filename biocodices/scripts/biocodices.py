#!/usr/bin/env python

import sys
import argparse
import yaml
from termcolor import colored
from os.path import expanduser, join, dirname

from biocodices import SequencerRun, Cohort, software_name
from biocodices.helpers.general import touch_all_the_logs, biocodices_logo


def call_variants_for_all_samples(base_dir):
    sequencer_run = SequencerRun(expanduser(base_dir))
    samples = sequencer_run.samples()
    cohort = Cohort(samples)

    print(biocodices_logo())
    print('Welcome to {}! Anlyzing reads for:'.format(software_name))
    print(colored(sequencer_run, 'green'))
    print('\nYou can follow the details of the process with:')
    # print('`tail -n0 -f {}/{{*/,}}*.log`\n'.format(sequencer_run.results_dir))
    print('`tail -n0 -f {}/*/*.log`\n'.format(sequencer_run.results_dir))

    dir_list = [sample.results_dir for sample in samples]
    touch_all_the_logs(sequencer_run.dir, dir_list)

    cohort.call_variants()

    print('\nDone! Bless your heart.\n')


if __name__ == '__main__':
    cli_yaml_path = join(dirname(__file__), 'cli_arguments.yml')
    with open(cli_yaml_path, 'r') as cli_yaml_file:
        cli = yaml.load(cli_yaml_file)

    parser = argparse.ArgumentParser(description=cli['description'])
    for flag, (long_name, action) in cli['args'].items():
        print(flag, long_name, action)

    #  parser = argparse.ArgumentParser(description=description)
    #  parser.add_argument('-j', '--joint-genotyping', action='store_true',
                        #  help=j_help)


    #  args = sys.argv
    #  if len(args) != 2:
        #  print('\nUsage:')
        #  print('{} <root dir of sequencer_run>'.format(__file__))
        #  sys.exit()

    #  base_dir = args[1]
    #  call_variants_for_all_samples(base_dir)
