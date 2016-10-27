from os import mkdir
from os.path import join, expanduser, abspath, basename, isdir
from glob import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from contextlib import redirect_stdout
from io import StringIO


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
        It will JSONify dicts and lists.
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

        columns_to_jsonify = []
        for column_name, series in df.iteritems():
            types = series.dropna().map(type).unique()
            if len(types) == 1 and types[0] in (list, dict):
                columns_to_jsonify.append(column_name)

        if columns_to_jsonify:
            # I need to copy the df to serialize some columns without modifying
            # the original dataframe that was passed as an argument.
            df = df.copy()
            for column_name in columns_to_jsonify:
                df[column_name] = df[column_name].map(json.dumps)

        df.to_csv(filepath, index=index, **kwargs)
        print('Written to', filepath)

        if columns_to_jsonify:
            print('JSONified columns: %s' % columns_to_jsonify)

        return filepath

    def read_csv(self, filename, subdir='results', **kwargs):
        """
        Read a CSV file in any of the Project's subdirectories
        (default='results'). It will try to parse as JSON the fields with
        dtype=np.object and convert them to Python objects.
        Extra arguments are passed to pd.read_csv().
        """
        filepath = join(self.dir, subdir, filename)
        print('Reading "{}"'.format(filename))
        df = pd.read_csv(filepath, **kwargs)

        for column_name, series in df.items():
            if not series.dtype == np.dtype('object'):
                # It's faster to skip non-object columns than to try and fail
                continue

            try:
                df[column_name] = series.fillna('""')\
                                        .map(json.loads)\
                                        .replace('', np.nan)
                print('  Parsed "{}" as JSON'.format(column_name))
            except (TypeError, json.JSONDecodeError):
                pass

        f = StringIO()
        with redirect_stdout(f):
            df.info(memory_usage='deep')

        info_lines = [line for line in f.getvalue().split('\n') if 'memory' in line]
        print('\n'.join(info_lines), '\n')
        return df

    def read_results_df(self, filename, **kwargs):
        return pd.read_csv(self.results_file(filename), **kwargs)

    def save_last_plot(self, filename):
        filepath = self.results_file(filename)
        plt.savefig(filepath, bbox_inches='tight')
        print('Written to', filepath)
