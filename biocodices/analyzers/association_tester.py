import pandas as pd
from os.path import join
from glob import glob

from biocodices.programs import Plink


class AssociationTester:
    def __init__(self, dataset):
        self.dataset = dataset

    def __repr__(self):
        return '<{} for dataset "{}">'.format(self.__class__.__name__,
                                              self.dataset.label)

    def read_plink_test(self, test_name, path):
        test_data = Plink.config('association_tests')[test_name]
        df = pd.read_table(path, sep='\s+')
        df.set_index(['CHR', 'SNP'], inplace=True)
        ret = pd.DataFrame({})
        if test_name == 'model_fisher':
            col = test_data['interesting_column']
            for model, df_model in df.groupby('TEST'):
                series = df_model[test_data['interesting_column']]
                series.name = '{}_{}_{}'.format(model, test_name, col)
                ret = ret.append(series)
            return ret.transpose()
        elif test_name in ['allelic_fisher', 'allelic_fisher_adjusted',
                           'allelic_chisq']:
            col = test_data['interesting_column']
            series = df[col]
            series.name = '{}_{}'.format(test_name, col)
            return series.to_frame()

    #  @staticmethod
    #  def _read_plink_model_table(path_to_table):
        #  tests = ['GENO', 'TREND', 'ALLELIC', 'DOM', 'REC']
        #  df = pd.read_table(path_to_table, sep='\s+')
        #  return df

    #  @staticmethod
    #  def _read_plink_assoc_fisher(path):
        #  df = pd.read_table(path)
        #  return df

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
