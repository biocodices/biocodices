import sklearn.decomposition
from math import sqrt
from pandas import DataFrame
from analyzers.base_pca import BasePCA
from helpers.helpers import percentage_fmt


class SklearnPCA(BasePCA):
    def __init__(self, dataset):
        self.dataset = dataset

    def run(self, normalize=True):
        """
        Receive a dataset with genotypes. The genotypes should be a DataFrame
        which index should at least have "sample" IDs (e.g. HG00096 ...) and
        ideally would also be a MultiIndex with "population" and
        "superpopulation" levels.

        Return:
        - a PCA object that responds to 'result' and 'explained_variance'
          * 'result' is a DataFrame with the same [Multi]Index but the
            eigenvalues as columns.
          * 'explained_variance' is a DataFrame with the explained variance
            per eigenvector.

        """
        genotypes = self.dataset.genotypes()
        if normalize:
            genotypes = genotypes.apply(self._normalize)

        # Leave only SNPs with genotype defined at every sample
        genotypes.dropna(axis=1, inplace=True)

        components_to_take = 15
        sklearn_pca = sklearn.decomposition.PCA()
        pca_result_matrix = sklearn_pca.fit_transform(genotypes.values)
        result = DataFrame(pca_result_matrix, index=genotypes.index)
        result = result.ix[:, :components_to_take-1]  # (Indexing in pandas)
        result.columns = ["PC{}".format(ix+1)
                          for ix in result.columns]
        result = result.sort_index()
        self.result = result

        variance = sklearn_pca.explained_variance_ratio_[:components_to_take]
        variance = DataFrame(variance, index=result.columns)
        variance.index.name = 'component'
        variance.columns = ['percentage']
        variance['percentage'] = variance['percentage'].map(percentage_fmt)
        self.explained_variance = variance
        self.write_result_csvs()  # Useful for d3 use

    def _normalize(self, series):
        # Taken from Patterson et al. 2006, doi:10.1371/journal.pgen.0020190
        mu = series.mean()
        p = mu/2
        q = 1 - p

        if mu == 0:
            return series - mu
        else:
            return (series - mu) / sqrt(p * q)
