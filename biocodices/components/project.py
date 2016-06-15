from os import makedirs
from os.path import join, expanduser, abspath, basename
from glob import glob


class Project:
    def __init__(self, base_dir):
        self.dir = join(abspath(expanduser(base_dir)))
        self.label = basename(self.dir)
        self._prepare_dirs()

    def __repr__(self):
        return '<{} "{}">'.format(self.__class__.__name__, self.label)

    def data_files(self, pattern=None):
        return [basename(fn) for fn in
                glob(join(self.data_dir, (pattern or '*')))]

    def result_files(self, pattern=None):
        return [basename(fn) for fn in
                glob(join(self.results_dir, (pattern or '*')))]

    def result_file(self, filename):
        return join(self.results_dir, filename)

    def data_file(self, filename):
        return join(self.data_dir, filename)

    def _prepare_dirs(self):
        self.results_dir = join(self.dir, 'results')
        self.data_dir = join(self.dir, 'data')
        makedirs(self.results_dir, exist_ok=True)
