from os import mkdir
from os.path import join, expanduser, abspath, basename, isdir
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt


class Project:
    def __init__(self, base_dir):
        self.dir = abspath(expanduser(base_dir))
        self.name = basename(self.dir)
        self.data_dir = join(self.dir, 'data')
        self.results_dir = join(self.dir, 'results')

        for directory in [self.dir, self.data_dir, self.results_dir]:
            if not isdir(directory):
                mkdir(directory)

    def __repr__(self):
        return '<{} "{}">'.format(self.__class__.__name__, self.name)

    def data_files(self, pattern=None):
        return [fn for fn in glob(join(self.data_dir, (pattern or '*')))]

    def results_files(self, pattern=None):
        return [fn for fn in glob(join(self.results_dir, (pattern or '*')))]

    def results_file(self, filename):
        return join(self.results_dir, filename)

    def data_file(self, filename):
        return join(self.data_dir, filename)

    def dump_df(self, df, filename, index=None, **kwargs):
        """
        Dump the dataframe to a CSV in the results dir with the given filename.
        Include '.tsv' in the filename to make it a TSV file!
        Extra **kwargs are passed to pandas.DataFrame.to_csv()
        """

        if 'sep' not in kwargs:
            if filename.endswith('.tsv'):
                kwargs['sep'] = '\t'
            else:
                kwargs['sep'] = ','

        if kwargs['sep'] == ',' and not filename.endswith('.csv'):
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
        df.to_csv(filepath, index=index, **kwargs)
        print('Written to', filepath)

        return filepath

    def read_results_df(self, filename, **kwargs):
        return pd.read_csv(self.results_file(filename), **kwargs)

    def save_last_plot(self, filename):
        filepath = self.results_file(filename)
        plt.savefig(filepath, bbox_inches='tight')
        print('Written to', filepath)
