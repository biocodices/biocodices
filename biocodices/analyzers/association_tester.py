import pandas as pd
from os.path import basename

from biocodices.programs import Plink
from biocodices.plotters import AssociationTestPlotter


class AssociationTester:
    def __init__(self, dataset):
        self.dataset = dataset
        self.plotter = AssociationTestPlotter(self)
        self.plink = Plink(self.dataset.label_path)
        self.tests = Plink.config('association_tests')

    def __repr__(self):
        return '<{} for dataset "{}">'.format(self.__class__.__name__,
                                              self.dataset.label)

    @staticmethod
    def available_tests():
        return list(Plink.config('association_tests').keys())

    def run(self, test_name, famfile=None):
        result_file = self.plink.association_test(test_name, famfile=famfile)
        result = self.read_plink_test(test_name, result_file)
        return result

    def read_plink_test(self, test_name, path):
        interesting_col = self.tests[test_name]['interesting_column']

        if test_name in ['five_models_fisher_for_case_control']:
            return self.read_model_file(path, interesting_col)

        if test_name in ['linear_regression_for_quanti_trait']:
            return self.read_qassoc_means_file(path, interesting_col)

    @staticmethod
    def read_model_file(path, interesting_col):
        df = pd.read_table(path, sep='\s+')
        df.set_index(['CHR', 'SNP'], inplace=True)
        ret = pd.DataFrame({})
        for model, df_model in df.groupby('TEST'):
            series = df_model[interesting_col]
            series.name = '{}_{}'.format(model, interesting_col)
            ret = ret.append(series)
        return ret.transpose()

    @staticmethod
    def read_qassoc_means_file(path, interesting_col):
        df = pd.read_table(path, sep='\s+')
        df.set_index(['CHR', 'SNP'], inplace=True)
        return df

        #  elif test_name in ['allelic_regression']:
            #  series = df[interesting_col]
            #  series.name = '{}_{}'.format(test_name, interesting_col)
            #  return series.to_frame()

    #  def plot_plink_results(self, results_file):
        #  df = self.read_plink_test(results_file)

    @staticmethod
    def dataframe_to_pheno_file(pheno_df, outfile):
        """
        Read phenotypes from a pandas dataframe (must have FID and IID columns)
        and write a .pheno file to use as PLINK input.
        (output tsv will have FID IID PHENO-1 PHENO-2 ... PHENO-N columns)
        """
        if not outfile.endswith('.pheno'):
            outfile += '.pheno'

        # Sort the columns for the .pheno file format
        sorted_df = pd.DataFrame({})
        sorted_df['FID'] = pheno_df['FID']
        sorted_df['IID'] = pheno_df['IID']
        rest_of_cols = list(set(pheno_df.columns) - set(['FID', 'IID']))
        sorted_df = pd.concat([sorted_df, pheno_df[rest_of_cols]], axis=1)
        sorted_df.sort_values(by=['FID', 'IID'], inplace=True)
        sorted_df.fillna('-9', inplace=True)

        sorted_df.to_csv(outfile, header=True, index=False)
        print('Written -> {}'.format(outfile))
        return outfile



    #  def run(self):
        #  self.result = {}
        #  tests = {}
        #  tests['UNADJ'] = self.dataset.plink.assoc(adjust=False)
        #  tests['ADJUST'] = self.dataset.plink.assoc(adjust=True)
        #  tests['model'] = self.dataset.plink.model()

        #  for test_name, test_filename in tests.items():
            #  df = pd.read_table(test_filename, sep='\s+').set_index('SNP')
            #  self.result[test_name] = df

        #  # Split the single table with all models into different tables
        #  for model_name, df in self.result['model'].groupby('TEST'):
            #  self.result[model_name] = df
        #  del(self.result['model'])

    # FIXME: needs rewriting
    #  def plot(self, ax=None, tests_to_plot=None):

        #  if ax is None:
            #  _, ax = plt.subplots(figsize=(15, 5))

        #  if tests_to_plot is None:
            #  tests_to_plot = self.available_tests

        #  self.plotter.draw_ax(ax, tests_to_plot)
        #  return ax
