from datetime import datetime
from os import makedirs
from os.path import join, expanduser, abspath
from glob import glob


class Project:
    def __init__(self, base_dir, label):
        # timestamp = datetime.now().strftime('%Y-%m-%d')
        # self.label = '{}__{}'.format(timestamp, label)
        self.label = label
        self.dir = join(abspath(expanduser(base_dir)), self.label)
        self._prepare_dirs()

    def __repr__(self):
        return '<{} "{}">'.format(self.__class__.__name__, self.label)

    @property
    def data_files(self, pattern=None):
        return glob(join(self.dir, (pattern or '*')))

    def _prepare_dirs(self):
        # Generate a plain text file with the date of the start of the project
        self.results_dir = join(self.dir, 'results')
        self.data_dir = join(self.dir, 'data')

        for directory in [self.dir, self.results_dir, self.data_dir]:
            makedirs(directory, exist_ok=True)


        with open(join(self.results_dir, 'project_date_of_creation'), 'w'):

