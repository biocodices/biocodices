import sys
import yaml
from functools import partial
from os import makedirs
from os.path import join, abspath, dirname, isfile
from itertools import product
from collections import OrderedDict
from datetime import datetime
from multiprocessing import Pool
from termcolor import colored

from biocodices import software_name
from biocodices.components.cohort import Cohort, EmptyCohortException
from biocodices.helpers import Stopwatch


class Pipeline:
    def __init__(self, log_filepath):
        self.tasks = OrderedDict()
        self.log_filepath = log_filepath

    def __repr__(self):
        return '<{} with tasks:\n{}>'.format(self.__class__.__name__,
                                             '\n'.join(self.tasks.keys()))

    def add_task(self, task, label):
        """Add tasks that should be either a method or a partial for later
        execution."""
        self.tasks[label] = task

    def add_multitask(self, task_group, task_group_label, n_processes=1):
        """task_group should be a dictionary of task_labels as keys and tasks
        methods or partials as values. The tasks will be run in the specified
        number of parallel processes."""
        if n_processes > 8:
            raise ValueError("I won't run more than 8 parallel processes.")

        task_group_label += ' (parallelized x {})'.format(n_processes)
        self.add_task(partial(self.run_task_group, task_group=task_group,
                              n_processes=n_processes), task_group_label)

    def run_task_group(self, task_group, n_processes=1):
        """This was meant as a pseudo-task that will actually run a group of
        tasks in parallel, so instead of queueing each task separatedly, this
        method will be queued in the pipeline by add_multitask() and this
        method will deal with the parallelization.
        """
        with Pool(processes=n_processes) as pool:
            for task_label, task in task_group.items():
                self.print_and_log_to_file(task_label)
                pool.apply_async(task)

            pool.close()  # Required call before pool.join()
            pool.join()  # Wait for all processes to finish before continuing

    def run(self):
        """Runs the whole pipeline, task by task, in a serial manner."""
        stopwatch = Stopwatch().start()
        options_in_effect = {k: v for k, v in vars(self.cli_args).items()
                             if v is not None}
        initial_message = 'Started the pipeline. Options in effect:\n\n' + \
            yaml.dump(options_in_effect, default_flow_style=False)
        self.print_and_log_to_file(initial_message, create=True)

        for task_label, task in self.tasks.items():
            self.print_and_log_to_file(task_label)
            task()

        self.total_time = stopwatch.stop()
        finish_message = 'Finished the pipeline. Total time: {}.'
        self.print_and_log_to_file(finish_message.format(self.total_time),
                                   print_it=False)

    def print_and_log_to_file(self, msg, create=False, print_it=True):
        open_mode = 'w' if create else 'a'
        with open(self.log_filepath, open_mode) as log:
            prefix = '[{}] '.format(self.timestamp())
            log.write(prefix + msg + '\n')
        if print_it:
            print(prefix + msg)

    @staticmethod
    def timestamp():
        return datetime.now().strftime('%H:%M:%S')


class PipelineCreator:
    def __init__(self, args):
        self.args = args

    def pre_pipeline(self):
        self._print_welcome_message()
        self.cohort = self._create_cohort(self.args.seq_dir)
        self._set_cohort_samples(keep=self.args.samples,
                                 skip=self.args.skip_samples)
        self._print_intro_information()
        self._touch_all_the_logs()

    def build_pipeline(self):
        """
        Builds a Pipeline storing the tasks to run in it. The pipeline will
        be stored in self until it's ran. It needs to have a cohort previously
        set.
        """
        pipeline = Pipeline(log_filepath=join(self.cohort.dir, 'pipeline.log'))
        pipeline.cli_args = self.args

        if self.args.trim_reads:
            task_group = OrderedDict()
            task_group_label = 'Trim and analyze reads'
            for sample in self.cohort.samples:
                task_label = sample.msg('Trim and analyze reads')
                task_group[task_label] = sample.analyze_and_trim_reads
            pipeline.add_multitask(task_group, task_group_label,
                                   n_processes=self.args.number_of_processes)
        if self.args.align_reads:
            for sample in self.cohort.samples:
                pipeline.add_task(sample.align_reads,
                                  sample.msg('Align reads'))
                pipeline.add_task(sample.process_alignment_files,
                                  sample.msg('Process alignment files'))
                pipeline.add_task(sample.alignment_metrics,
                                  sample.msg('Alignment metrics'))

        if self.args.create_vcfs:
            task_group = OrderedDict()
            task_group_label = 'Haplotype Caller'
            for sample in self.cohort.samples:
                task_label = sample.msg('Create gVCF')
                task_group[task_label] = sample.create_gvcf
            pipeline.add_multitask(task_group, task_group_label,
                                   n_processes=self.args.number_of_processes)

        if self.args.plot_metrics:
            pipeline.add_task(self.cohort.plot_alignment_metrics,
                              self.cohort.msg('Plot alignment metrics'))
            pipeline.add_task(self.cohort.median_coverages,
                              self.cohort.msg('Compute median coverages'))

        if self.args.joint_genotyping:
            pipeline.add_task(self.cohort.joint_genotyping,
                              self.cohort.msg('Joint genotyping'))
        if self.args.hard_filtering:
            pipeline.add_task(partial(self.cohort.apply_filters_to_vcf,
                                      self.cohort.unfiltered_vcf),
                              self.cohort.msg('Hard filtering'))
            for sample in self.cohort.samples:
                pipeline.add_task(
                    partial(self.cohort.subset_samples,
                            self.cohort.filtered_vcf,
                            [sample.id],
                            sample.filtered_vcf),
                    self.cohort.msg('Subset from multisample VCF'))

        self.pipeline = pipeline
        return pipeline

        # self.printlog('Plot some alignment metrics for the cohort.')
        # self.printlog('Compute median coverage of the cohort.')
        # self.printlog('Joint genotyping.')
        # self.printlog('Hard filtering the multisample VCF')
        # self.printlog('Split the multisample VCF per sample')

    def post_pipeline(self):
        self._print_exit_message(total_time=self.pipeline.total_time)

    @staticmethod
    def _create_cohort(seq_dir):
        try:
            return Cohort(seq_dir)
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
        print(cls._biocodices_logo())
        print('Welcome to {}! Reading the data directory...'.format(software_name))

    @staticmethod
    def _biocodices_logo():
        path = abspath(join(dirname(__file__), '../scripts/logo.txt'))
        with open(path, 'r') as logo_file:
            logo_string = logo_file.read()
        return logo_string

    def _print_intro_information(self):
        print(colored(self.cohort, 'green'))
        print('\nYou can follow the details of the process with:')
        # other option: `tail -n0 -f {}/{{*/,}}*.log`
        print('`tail -n0 -f {}/{{*/,}}*.log`\n'.format(self.cohort.results_dir))

    def _touch_all_the_logs(self):
        # I wrote this just to be able to run a `tail -f *.log` on every log
        # during the variant calling, even for logs that don't yet exist but
        # that would later be created, so what I do is just creating them
        # beforehand. It's a necessarily hardcoded list, I guess:
        samples_dirs = [sample.results_dir for sample in self.cohort.samples]
        for d in samples_dirs:
            makedirs(d, exist_ok=True)

        lognames_file = abspath(join(dirname(__file__),
                                     '../scripts/log_filenames.txt'))
        with open(lognames_file, 'r') as f:
            log_filenames = [l.strip() for l in f.readlines()]
        log_filepaths = [abspath(join(d, fn + '.log'))
                         for d, fn in product(samples_dirs, log_filenames)]
        for log_filepath in log_filepaths:
            if not isfile(log_filepath):
                open(log_filepath, 'a').close()

    @staticmethod
    def _print_exit_message(total_time):
        print('\nThe whole process took {}.'.format(total_time))
        print('Done! Bless your heart.\n')
