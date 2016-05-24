import pandas as pd

from os.path import basename

from biocodices.programs.plink import Plink
from biocodices.helpers.config import Config


class Panel:
    def __init__(self, bim_filepath):
        """
        Reads SNPs info from a bim file. Assumes a file format like
        '<panel_name>.<samplegroup_name>.bim'.
        TODO: Reads extra info for those SNPs in a database.
        """
        self.file = bim_filepath
        self.label = basename(self.file.split('.')[0])
        # self.name = Config('names').get(self.label, self.label)
        self.name = self.label
        self.markers = self._read_bim(self.file)
        self.snps = self.markers  # Legacy alias, TODO: remove it

        # self.snps_file = self.info_path + '.snps'

        # info_file = self.info_path + '.csv'
        #  if isfile(info_file):
            #  self.extra_info = self.read_info(info_file)

        self.rs_ids = self.snps.index.values  # Redundant, but handy shortcut
        self.name = self._generate_name()

    def __repr__(self):
        return '<Panel {}>'.format(self.name)

    def __len__(self):
        return len(self.snps)

    def _generate_name(self):
        return '{0} · {1:,} SNPs'.format(self.label, len(self))

    @staticmethod
    def _read_bim(filename):
        return pd.read_table(filename, names=Plink.bim_fields,
                             index_col="rs_id")

    #  @staticmethod
    #  def read_info(filename):
        #  return pd.read_csv(filename, index_col="rs_id")

    #  @classmethod
    #  def available_panels(cls, source_label=None):
        #  bim_files = glob(join(cls.base_dir(), '*.bim'))
        #  if source_label is not None:
            #  glob_expr = join(Source(source_label).panels_dir, '*.bim')
            #  bim_files = glob(glob_expr)
        #  panel_labels = [basename(f).replace('.bim', '') for f in bim_files]
        #  return sorted(panel_labels)
