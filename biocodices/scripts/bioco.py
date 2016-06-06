#!/usr/bin/env python

import argparse
import yaml
from os.path import join, dirname, expanduser, abspath
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
# This prevents matplotlib raising an exception when running biocodices on a
# remote server with no X. This line has to be executed before importing
# pyplot, so we have to run it before any biocodices code is imported.

from biocodices.components import PipelineCreator


if __name__ == '__main__':
    cli_info_path = join(dirname(__file__), 'cli_arguments.yml')
    with open(cli_info_path, 'r') as cli_info_file:
        cli = yaml.load(cli_info_file)

    parser = argparse.ArgumentParser(description=cli['description'])
    for flag, options in cli['args'].items():
        long_name, action, required, default = options
        parser.add_argument('-' + flag, '--' + long_name, action=action,
                            required=required, help=cli['help_texts'][flag],
                            default=default)

    args = parser.parse_args()
    args.seq_dir = abspath(expanduser(args.seq_dir))
    args.number_of_processes = int(args.number_of_processes)
    if args.complete_pipeline:
        for attr in ['trim_reads', 'align_reads', 'create_vcfs',
                     'plot_metrics', 'joint_genotyping', 'hard_filtering']:
            setattr(args, attr, True)

    pipeline_creator = PipelineCreator(args)
    pipeline_creator.pre_pipeline()
    pipeline = pipeline_creator.build_pipeline()
    pipeline.run()
    pipeline_creator.post_pipeline()
