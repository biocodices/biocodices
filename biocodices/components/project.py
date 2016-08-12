from os import mkdir
from os.path import join, expanduser, abspath, basename, isdir
from glob import glob
import pandas as pd


class Project:
    def __init__(self, base_dir):
        self.dir = abspath(expanduser(base_dir))
        self.id = basename(self.dir)
        self.data_dir = join(self.dir, 'data')
        self.results_dir = join(self.dir, 'results')

        for directory in [self.dir, self.data_dir, self.results_dir]:
            if not isdir(directory):
                mkdir(directory)

    def __repr__(self):
        return '<{} "{}">'.format(self.__class__.__name__, self.id)

    def data_files(self, pattern=None):
        return [fn for fn in glob(join(self.data_dir, (pattern or '*')))]

    def results_files(self, pattern=None):
        return [fn for fn in glob(join(self.results_dir, (pattern or '*')))]

    def results_file(self, filename):
        return join(self.results_dir, filename)

    def data_file(self, filename):
        return join(self.data_dir, filename)

    def dump_df(self, df, filename, index=None):
        """
        Dump the dataframe to a CSV in the results dir with the given filename.
        """
        if not filename.endswith('.csv'):
            filename += '.csv'

        if index is None:
            # Don't include the index if it's just ordered numbers
            # This guessing can be overriden by specifying 'index'
            num_index_types = [
                pd.indexes.range.RangeIndex,
                pd.indexes.numeric.Int64Index
            ]
            index = (type(df.index) not in num_index_types)

        filepath = self.results_file(filename)
        df.to_csv(filepath, index=index)
        print('Written to', filepath)

        return filepath
