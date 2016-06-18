import yaml
from functools import partial
from multiprocessing import Pool
from collections import OrderedDict

from biocodices.helpers import Stopwatch, timestamp


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
                # self.print_and_log_to_file(task_label)
                pool.apply_async(task)

            pool.close()  # Required call before pool.join()
            pool.join()  # Wait for all processes to finish before continuing

    def run(self):
        """Runs the whole pipeline, task by task, in a serial manner."""
        stopwatch = Stopwatch().start()
        options_in_effect = {k: v for k, v in self.cli_args.items() if v}
        initial_message = 'Started the pipeline. Options in effect:\n\n' + \
            yaml.dump(options_in_effect, default_flow_style=False)
        self.print_and_log_to_file(initial_message, create=True)

        for task_label, task in self.tasks.items():
            self.print_and_log_to_file(task_label)
            result = task()
            print('result:', result)

        self.total_time = stopwatch.stop()
        finish_message = 'Finished the pipeline. Total time: {}.'
        self.print_and_log_to_file(finish_message.format(self.total_time))

    def print_and_log_to_file(self, msg, create=False, print_it=True):
        open_mode = 'w' if create else 'a'
        with open(self.log_filepath, open_mode) as log:
            prefix = '[{}] '.format(timestamp())
            log.write(prefix + msg + '\n')
        if print_it:
            print(prefix + msg)
