import pandas as pd
from os.path import expanduser, isfile, basename, dirname

from biocodices.analyzers import (SmartPCA, SklearnPCA, Admixture,
                                  QualityControl)
from .panel import Panel
from .sample_group import SampleGroup
from biocodices.programs.plink import Plink


class Dataset:
    def __init__(self, plink_label_path):
        self.label_path = expanduser(plink_label_path)
        self.dir = dirname(self.label_path)
        self._check_plinkfiles()
        self.label = basename(self.label_path)

    def __repr__(self):
        tmpl = '[ Dataset from:\n  {}\n  {} ]'
        return tmpl.format(self.panel, self.samplegroup)

    @property
    def samplegroup(self):
        return SampleGroup(self.fam)

    @property
    def panel(self):
        return Panel(self.bim)

    @property
    def markers(self):
        return self.panel.markers

    @property
    def samples(self):
        return self.samplegroup.samples

    @property
    def genotypes(self):
        # if not isfile(self.traw):
        # ^ I re-create the .traw on demand now because i ran into situations
        # where the dataset got 'stuck' with an old traw file that was result
        # of old runs of the code in the same directory, so while the
        # .genotypes property showed something, the .bed actually had
        # something else.

        Plink(self.label_path).make_traw()  # Create just to read genotypes from
        return self._read_traw(self.traw)

    def quality_control(self):
        """
        Run a battery of quality control tests and store the output in a
        QualityControl object that responds to #result and #plot().
        """
        quali = QualityControl(dataset=self)
        quali.run()
        return quali

    def apply_pre_test_filters(self, mind=0.1, maf=0.05, geno=0.1, hwe=0.0001,
                               out=None):
        """
        Apply some PLINK filters before running the association tests:
            mind: leave out samples with a call rate below the threshold.
            geno: leave out variants with a call rate exceding (?) the
            threshold.
            maf: leave out variants with a minor allele frequency below the
            threshold.
            hew: leave out variants which have a Hardy-Weiberg equilibrium
            exact test p-value below the threshold.
        Returns a new (filtered) dataset.
        """
        out = out or (self.label_path + '_filtered_for_tests')
        Plink(self.label_path).pre_tests_filter(mind, maf, geno, hwe, out)
        return Dataset(out)

    def pca(self, implementation='smartpca', overwrite=False,
            normalize=True, args={}):
        """
        Computes a Principal Components Analysis with the genotypes in this
        dataset. Returns a PCA object that responds to #result and #plot().
        """
        # SmartPCA needs a ped file to read the genotypes
        # Re-create it each time to make sure it's updated with any changes
        # to the famfile in this dataset.
        Plink(self.label_path).make_ped()

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
        return '{}.{}'.format(self.label_path, extension)

    def _check_plinkfiles(self):
        self.bed = self._filename('bed')
        self.ped = self._filename('ped')
        self.bim = self._filename('bim')
        self.fam = self._filename('fam')
        self.traw = self._filename('traw')

        if not isfile(self.bed):
            if not isfile(self.ped):
                raise FileNotFoundError('No bed or ped found there.')
            Plink.make_bed_from_ped(self.label_path)

    def _read_traw(self, filepath):
        df = pd.read_table(filepath, index_col='SNP')
        # ^ df is a genotypes matrix with values 0, 1, 2 for alt allele dosage.
        df.drop(['CHR', '(C)M', 'POS', 'COUNTED', 'ALT'], axis=1, inplace=True)
        df = df.transpose()
        df.columns.name = 'rs_id'
        fids = [iid_fid.split('_')[0] for iid_fid in df.index]
        iids = [iid_fid.split('_')[1] for iid_fid in df.index]
        df['IID'], df['FID'] = iids, fids
        df.reset_index(drop=True, inplace=True)
        df.set_index(['FID', 'IID'], inplace=True)
        # merge the genotype data (.bed) with the samples data (.fam)
        df = self.samplegroup.samples.join(df)
        return df.sort_index()
