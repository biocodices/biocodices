import pandas as pd
from os.path import expanduser, isfile, basename, dirname

from biocodices.analyzers.smart_pca import SmartPCA
from biocodices.analyzers.sklearn_pca import SklearnPCA
from biocodices.analyzers.admixture import Admixture
from biocodices.analyzers.association_test import AssociationTest
from biocodices.analyzers.quality_control import QualityControl
from .panel import Panel
from .sample_group import SampleGroup
from biocodices.helpers.plink import Plink


class Dataset:
    def __init__(self, plink_label_path):
        self.path_label = expanduser(plink_label_path)
        self.dir = dirname(self.path_label)
        self._check_plinkfiles()

        self.label = basename(self.path_label)
        self.panel = Panel(self.bimfile)
        self.samplegroup = SampleGroup(self.famfile)
        self.genotypes = self._read_traw(self.trawfile)

    def __repr__(self):
        tmpl = '[ Dataset from:\n  {}\n  {} ]'
        return tmpl.format(self.panel, self.samplegroup)

    def genotypes(self):
        raise Exception('Not yet implemented')

    def quality_control(self):
        """
        Run a battery of quality control tests and store the output in a
        QualityControl object that responds to #result and #plot().
        """
        quali = QualityControl(dataset=self)
        quali.run()
        return quali

    def association_tests(self):
        """
        Run a battery of association tests per SNP and store the output in
        an AssociationTest object that responds to #result and #plot().
        """
        assoc = AssociationTest(dataset=self)
        assoc.run()
        return assoc

    def pca(self, implementation='smartpca', overwrite=False,
            normalize=True, args={}):
        """
        Computes a Principal Components Analysis with the genotypes in this
        dataset. Returns a PCA object that responds to #result and #plot().
        """
        if implementation == 'smartpca':
            pca = SmartPCA(dataset=self)
            pca.run(overwrite=overwrite, args=args)
            return pca
        elif implementation == 'sklearn':
            pca = SklearnPCA(self)
            pca.run(normalize=normalize)
            return pca

    def admixture(self, Ks, cores=4, infer_components=True, overwrite=False):
        """
        Run admixture for a list of K values. You may specify the number of
        CPU cores for admixture to use. Set 'overwrite=True' to make it ignore
        existing results files if it finds them; otherwise, it will read them
        instead of rerunning the analysis, to save time.
        Returns an admixture object that responds to #result and #plot().
        """
        admixture = Admixture(self)
        admixture.run(Ks, cores, infer_components, overwrite)
        return admixture

    def _filename(self, extension):
        return '{}.{}'.format(self.path_label, extension)

    def _check_plinkfiles(self):
        self.bedfile = self._filename('bed')
        self.pedfile = self._filename('ped')
        self.bimfile = self._filename('bim')
        self.famfile = self._filename('fam')
        self.trawfile = self._filename('traw')

        if not isfile(self.bedfile):
            if not isfile(self.pedfile):
                raise FileNotFoundError('No befile or pedfile found there.')
            Plink.make_bed_from_ped(self.path_label)

        self.plink = Plink(self.path_label)
        if not isfile(self.pedfile):
            self.plink.make_ped()  # Created for SmartPCA

        if not isfile(self.trawfile):
            self.plink.make_traw()  # Created to read genotypes from

    def _read_traw(self, filepath):
        df = pd.read_table(filepath, index_col='SNP')
        df.drop(['CHR', '(C)M', 'POS', 'COUNTED', 'ALT'], axis=1, inplace=True)
        # fids = [iid_fid.split('_')[0] for iid_fid in df.columns]
        iids = [iid_fid.split('_')[1] for iid_fid in df.columns]
        df.columns = iids  # Use individual IDs as columns!
        df.columns.name, df.index.name = 'sample', 'rs_id'
        df = self.samplegroup.samples.join(df.transpose())
        multi_index = ['sample']
        df = df.reset_index().set_index(multi_index)
        return df.sort_index()
