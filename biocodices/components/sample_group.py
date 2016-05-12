import pandas as pd

# from helpers.config import Config
from biocodices.helpers.plink import Plink


class SampleGroup:
    def __init__(self, fam_filepath):
        """
        Reads samples info from a fam file. Assumes a file format like
        '<panel_name>.<samplegroup_name>.fam'.
        TODO: Gather info for each sample in a database.
        """
        self.file = fam_filepath
        self.label = self.file.split('.')[1]
        # self.name = Config('names').get(self.label, self.label)
        self.samples = self._read_fam()

        # handy shortcuts
        self.sample_ids = self.samples.index.values
        # self.populations = self.samples['population'].unique()
        # self.regions = self.samples['region'].unique()

    def summary(self):
        n_families = len(self.samples['family'].unique())
        male_mask = self.samples['sexcode'] == 1
        n_males = len(self.samples[male_mask])
        n_females = len(self.samples[~male_mask])
        n_phenotypes = len(self.samples['phenotype'].unique())

        tmpl =  '{} samples, {} families, {} males & {} females, {} phenotypes'
        return tmpl.format(len(self.samples), n_families, n_males, n_females,
                           n_phenotypes)

    def __len__(self):
        return len(self.samples)

    def __repr__(self):
        tmpl = '<SampleGroup {} · {} samples>'
        return tmpl.format(self.label, len(self))

    def _read_fam(self):
        return pd.read_table(self.file, names=Plink.fam_fields(),
                             index_col="sample", sep='\s+')

    #  def _write_clusters_files(self, level='population'):
        #  """
        #  Write a file with three columns: FID, IID, and cluster (pop or region).
        #  FID and IID are usually the same, the clusters field is what matters
        #  for later use with plink --fst.
        #  """
        #  clusters = self.samples.reset_index()
        #  clusters["FID"] = clusters["sample"]
        #  # FIXME: FID=IID might not be ok if there's family info in the samples
        #  clusters = clusters[["FID", "sample", level]]
        #  filename = "{}.{}.clusters".format(self.label, level)
        #  filepath = join(self.base_dir, filename)
        #  clusters.to_csv(filepath, sep="\t", header=None, index=False)

        #  return filepath

    #  def clusters_filepath(self, level):
        #  filename = '{}.{}.clusters'.format(self.label, level)
        #  return join(self.base_dir, filename)
