import pandas as pd
from os.path import basename

from biocodices.programs import Plink
from biocodices.plotters import AssociationTestPlotter


class AssociationTester:
    def __init__(self, dataset):
        self.dataset = dataset
        self.plotter = AssociationTestPlotter(self)
        self.plink = Plink(self.dataset.label_path)
        self.result_files = {}
        self.results = {}

    def __repr__(self):
        return '<{} for dataset "{}">'.format(self.__class__.__name__,
                                              self.dataset.label)
    @property
    def tests(self):
        return Plink.config('association_tests')

    @staticmethod
    def available_tests():
        return list(Plink.config('association_tests').keys())

    def run(self, test_name, famfile=None):
        result_file = self.plink.association_test(test_name, famfile=famfile)
        result = self.read_plink_test(test_name, result_file)
        return result

    def read_plink_test(self, test_name, path):
        interesting_col = self.tests[test_name]['interesting_column']
        test_label = basename(path)

        if test_name in ['linear_regression_for_quanti_trait',
                         'trend_model_with_permutations',
                         'dom_model_with_permutations',
                         'rec_model_with_permutations',
                         'gen_model_with_permutations',
                         'allelic_model_with_permutations']:
            p_values_df = self.read_plink_results(path, interesting_col)
            self.result_files[test_label] = p_values_df
            self.results[test_label] = p_values_df
            return p_values_df

        # These are special cases since for each SNP tested there will be
        # five rows, one for each model tested.
        if test_name in ['five_models_fisher_for_case_control',
                         'five_models_chisq_for_case_control']:
            self.result_files[test_label] = path
            all_models, p_values_df = self.read_model_file(path,
                                                           interesting_col)
            self.results[test_label] = all_models
            return p_values_df

        msg = "I don't know how to read results for the test '{}'."
        raise ValueError(msg.format(test_name))

    @staticmethod
    def read_plink_results(path, interesting_col):
        df = pd.read_table(path, sep='\s+')
        df.set_index(['CHR', 'SNP'], inplace=True)
        return df

    @classmethod
    def read_model_file(cls, path, interesting_col):
        df = pd.read_table(path, sep='\s+')
        df.set_index(['CHR', 'SNP'], inplace=True)
        cls._add_cases_controls_columns(df)
        ret = pd.DataFrame({})
        for model, df_model in df.groupby('TEST'):
            series = df_model[interesting_col]
            series.name = '{}_{}'.format(model, interesting_col)
            ret = ret.append(series)
        return df, ret.transpose()

    @staticmethod
    def _add_cases_controls_columns(df):
        def sum_genotypes(genotype_cell):
            # genotype_cell is a string like 6/12/8
            return sum([int(n) for n in genotype_cell.split('/')])

        df.loc[:, 'n_cases'] = df['AFF'].map(sum_genotypes)
        df.loc[:, 'n_controls'] = df['UNAFF'].map(sum_genotypes)
        df.loc[:, 'cc_ratio'] = df['n_cases'] / df['n_controls']

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
