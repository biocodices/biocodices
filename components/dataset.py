import pandas as pd
from os.path import expanduser, isfile, basename

from components.panel import Panel
from components.sample_group import SampleGroup
from helpers.plink import Plink


class Dataset:
    def __init__(self, plink_label_path):
        self.path_label = expanduser(plink_label_path)
        self.label = basename(self.path_label)
        self._check_plinkfiles()
        self.plink = Plink(self.path_label)
        self.panel = Panel(self._filename('bim'))
        self.samplegroup = SampleGroup(self._filename('fam'))

    def __repr__(self):
        tmpl = '[ Dataset from:\n  {}\n  {} ]'
        return tmpl.format(self.panel, self.samplegroup)

    def assoc(self):
        self.plink.assoc(adjust=False)
        return self._read_assoc_output()

    def _read_assoc_output(self):
        filepath = self._filename('assoc')
        df = pd.read_table(filepath)
        return df

    def _filename(self, extension):
        return '{}.{}'.format(self.path_label, extension)

    def _check_plinkfiles(self):
        if isfile(self._filename('bed')):
            return

        if isfile(self._filename('ped')):
            Plink.make_bed_from_ped(self.path_label)
            return

        raise FileNotFoundError('No befile or pedfile found there.')
