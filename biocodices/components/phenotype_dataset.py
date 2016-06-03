from os.path import join
import pandas as pd
import seaborn as sns


class PhenotypeDataset:
    MISSING = 'NA'  # Plink takes non-numeric cells as a missing values
    BINARY_MISSING = -9  # Plink's code for missing value among case-controls

    def __init__(self, phenos_dataframe):
        """
        Initialize the class with a pandas DataFrame with:
            - a [multi]index with the ids of the samples
            - one column per phenotype
        Phenotypes can be binary, quantitative or both.
        """
        self.raw = phenos_dataframe
        self.all, self.binary, self.quanti = self._parse(self.raw)
        self.cases_vs_controls = self._cases_vs_controls(self.binary)
        self.missing_quanti = self._missing_quanti(self.quanti)

    def plot_cases_vs_controls(self, **kwargs):
        """
        Plots a stacked barchart of cases vs controls vs missing values. Any
        keyword arguments will be passed directly to pandas' plot function.
        """
        sns.set_style('white')
        ax = self.cases_vs_controls.plot(
            kind='barh', stacked=True, color=sns.color_palette(), **kwargs)
        ax.axvline(len(self.binary)/2, linestyle='dashed', color='indigo')
        ax.legend(bbox_to_anchor=(1.5, 1))
        sns.despine(left=True)
        return ax

    def plot_quanti_missing(self, **kwargs):
        """
        Plots a stacked barchart of cases vs controls vs missing values. Any
        keyword arguments will be passed directly to pandas' plot function.
        """
        sns.set_style('white')
        ax = self.missing_quanti.plot(
            kind='barh', stacked=True, color=sns.color_palette(), **kwargs)
        ax.legend(bbox_to_anchor=(1.5, 1))
        ax.axvline(len(self.quanti)/2, linestyle='dashed', color='indigo')
        sns.despine(left=True)
        return ax

    def write_a_fam_file_per_phenotype(self, dest_dir, original_fam_df):
        new_famfiles = {}
        for phenotype_name, phenotypes in self.all.iteritems():
            new_fam_filepath = '{}.fam'.format(phenotype_name)
            new_fam_filepath = join(dest_dir, new_fam_filepath)
            new_fam_df = original_fam_df
            # ^ Take the samples df from the datasets object;
            # it has the same columns than plink's fam files.

            if phenotype_name in self.binary.columns:
                # Convert a binary phenotype (0 control, 1 case, -9 missing)
                # to plink's fam file format (1 control, 2 case, -9 missing)
                new_fam_df['phenotype'] = phenotypes.map({0: 1, 1: 2, -9: -9})
            elif phenotype_name in self.quanti.columns:
                new_fam_df['phenotype'] = phenotypes
            else:
                msg = 'This pheno is neither binary nor quanti? {}'
                print(msg.format(phenotype_name))

            new_fam_df.to_csv(new_fam_filepath, sep='\t', header=False)
            new_famfiles[phenotype_name] = new_fam_filepath

        return new_famfiles

    def _parse(self, df):
        all_phenos = df.apply(self._int_if_binary)
        binary_columns = [col for col, series in df.items()
                          if self._is_this_series_binary(series)]
        binary = all_phenos[binary_columns]
        quanti = all_phenos.drop(binary_columns, axis=1)
        return (all_phenos, binary, quanti)

    @staticmethod
    def _is_this_series_binary(series):
        non_nan_values = set(series.dropna().unique())
        binary = len(non_nan_values) < 3
        return binary

    @classmethod
    def _int_if_binary(cls, series):
        # Transform to integers and -9 (missing) if the phenotype is binary
        # Leave floats and 'NA' (missing) if the phenotytpe is continuous

        if cls._is_this_series_binary(series):
            return series.fillna(cls.BINARY_MISSING).map(int)
        else:
            return series.map(float).fillna(cls.MISSING)

    @staticmethod
    def _cases_vs_controls(binary_phenos):
        cases_vs_controls = pd.DataFrame({})
        for phenotype, series in binary_phenos.items():
            count = series.value_counts()
            count.name = phenotype
            cases_vs_controls = cases_vs_controls.append(count)
        cases_vs_controls.sort_values(cases_vs_controls.columns[0], inplace=True)
        cases_vs_controls.fillna(0, inplace=True)
        return cases_vs_controls

    @classmethod
    def _missing_quanti(cls, quanti_phenos):
        def count_missing(series):
            return series[series == cls.MISSING].count()

        def count_present(series):
            return series[series != cls.MISSING].count()

        counts = pd.DataFrame({'missing': quanti_phenos.apply(count_missing),
                               'present': quanti_phenos.apply(count_present)})
        counts.sort_values(by=counts.columns[0], inplace=True)
        return counts
