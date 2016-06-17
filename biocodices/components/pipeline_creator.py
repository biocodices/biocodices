import os
import sys
from functools import partial
from collections import OrderedDict
from termcolor import colored

from biocodices import software_name
from biocodices.helpers import all_log_filepaths, logo, timestamp
from biocodices.components.cohort import Cohort, EmptyCohortException
from biocodices.components import Pipeline


class PipelineCreator:
    def __init__(self, arguments):
        """
        Arguments to be read from command line. Arguments must include a
        base directory to instatiate a Cohort and look for samples.
        """
        self.args = arguments

    def pre_pipeline(self):
        self._print_welcome_message()
        self.cohort = self._create_cohort(self.args['--base-dir'])
        self._set_cohort_samples(keep=self.args['--samples'],
                                 skip=self.args['--skip-samples'])
        self._print_intro_information()
        self._touch_all_the_logs()

    def build_pipeline(self):
        """
        Builds a Pipeline storing the tasks that will run.
        Once the setup is complete, the pipeline is run.
        """
        pipeline_filename = 'pipeline_{}.log'.format(timestamp(sep='',
                                                               date=True))
        pipeline = Pipeline(log_filepath=os.path.join(self.cohort.dir,
                                                      pipeline_filename))
        pipeline.cli_args = self.args
        n_processes = self.args['--parallel']

        if self.args['--trim-reads']:
            task_group = OrderedDict()
            for sample in self.cohort.samples:
                task_label = sample.msg('Trim and analyze reads')
                task_group[task_label] = sample.analyze_and_trim_reads

            pipeline.add_multitask(task_group,
                                   self.cohort.msg('Trim and analyze reads'),
                                   n_processes=n_processes)

        if self.args['--align-reads']:
            # I don't parallelize this step since the alignment already
            # works in parallel. Plus it takes a humongous amount of RAM.
            task_group = OrderedDict()
            for sample in self.cohort.samples:
                task_label = sample.msg('Align reads')
                task_group[task_label] = sample.align_reads

            # Limit the parallelization in this step to a max of 2 processes
            # since the aligner is a memory hog.
            pipeline.add_multitask(task_group,
                                   self.cohort.msg('Align reads'),
                                   n_processes=min(n_processes, 2))

            task_group = OrderedDict()
            for sample in self.cohort.samples:
                task_label = sample.msg('Process alignment files')
                task_group[task_label] = sample.process_alignment_files
            pipeline.add_multitask(task_group,
                                   self.cohort.msg('Alignment processing'),
                                   n_processes=n_processes)

            task_group = OrderedDict()
            for sample in self.cohort.samples:
                task_label = sample.msg('Alignment metrics')
                task_group[task_label] = sample.alignment_metrics
            pipeline.add_multitask(task_group,
                                   self.cohort.msg('Alignment metrics'),
                                   n_processes=n_processes)

        if self.args['--create-vcf']:
            task_group = OrderedDict()
            for sample in self.cohort.samples:
                task_label = sample.msg('Create gVCF')
                task_group[task_label] = sample.create_gvcf
            pipeline.add_multitask(task_group,
                                   self.cohort.msg('Haplotype Caller'),
                                   n_processes=n_processes)

        if self.args['--metrics']:
            pipeline.add_task(self.cohort.plot_alignment_metrics,
                              self.cohort.msg('Plot alignment metrics'))
            pipeline.add_task(self.cohort.median_coverages,
                              self.cohort.msg('Compute median coverages'))

        if self.args['--joint-genotyping']:
            pipeline.add_task(self.cohort.joint_genotyping,
                              self.cohort.msg('Joint genotyping'))

        if self.args['--hard-filtering']:
            pipeline.add_task(partial(self.cohort.apply_filters_to_vcf,
                                      self.cohort.unfiltered_vcf),
                              self.cohort.msg('Hard filtering'))

            task_group = OrderedDict()
            for sample in self.cohort.samples:
                task_label = sample.msg('Subset from multisample VCF')
                task = partial(self.cohort.subset_samples,
                               self.cohort.filtered_vcf,
                               [sample.id], sample.filtered_vcf)
                task_group[task_label] = task
            pipeline.add_multitask(task_group,
                                   self.cohort.msg('Subset from multisample VCF'),
                                   n_processes=n_processes)

        #  if self.args['--annotation']:


        self.pipeline = pipeline
        return pipeline

    def post_pipeline(self):
        self._print_exit_message(total_time=self.pipeline.total_time)

    @staticmethod
    def _create_cohort(base_dir):
        try:
            return Cohort(base_dir)
        except EmptyCohortException as error:
            print(error, '\n')
            sys.exit()

    def _set_cohort_samples(self, keep, skip):
        # Keep all by default
        keep_ids = (keep and keep.split(',')) or \
                   [s.id for s in self.cohort.samples]
        # Skip none by default
        skip_ids = (skip and skip.split(',')) or []

        for sample_id in keep_ids + skip_ids:
            if sample_id not in [sample.id for sample in self.cohort.samples]:
                msg = '{} not found in this cohort.'
                print(msg.format(sample_id))
                sys.exit()

        self.cohort.samples = [sample for sample in self.cohort.samples
                               if sample.id in keep_ids and
                               sample.id not in skip_ids]

    @classmethod
    def _print_welcome_message(cls):
        print(logo())
        print('Welcome to {}! Reading the data directory...'.format(software_name))

    def _print_intro_information(self):
        print(colored(self.cohort, 'green'))
        print('\nYou can follow the details of the process with:')
        print('`tail -n0 -f {}/{{*/,}}*.log`\n'.format(self.cohort.results_dir))

    def _touch_all_the_logs(self):
        # I wrote this just to be able to run a `tail -f *.log` on every log
        # during the variant calling, even for logs that don't yet exist but
        # that would later be created. What I do is just creating the empty
        # logs beforehand. It's a necessarily hardcoded list, I guess:
        samples_dirs = [sample.results_dir for sample in self.cohort.samples]
        for d in samples_dirs:
            os.makedirs(d, exist_ok=True)

        for log_filepath in all_log_filepaths(base_dir=self.cohort.results_dir,
                                              samples_dirs=samples_dirs):
            if not os.path.isfile(log_filepath):
                open(log_filepath, 'a').close()

    @staticmethod
    def _print_exit_message(total_time):
        print('\nThe whole process took {}.'.format(total_time))
        print('Done! Bless your heart.\n')
