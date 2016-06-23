import re
from glob import glob
from os import makedirs
from os.path import join, basename, dirname
from shutil import move
import pandas as pd

from biocodices.components import Sample, Project
from biocodices.variant_calling import VcfMunger
from biocodices.helpers import plural
from biocodices.plotters import AlignmentMetricsPlotter


class Cohort(Project):
    def __init__(self, base_dir):
        """
        Expects the path to a directory that will have a 'data' subdirectory
        with fastq forward (R1) and reverse (R2) files from samples in it.
        """
        super(self.__class__, self).__init__(base_dir)
        self._move_sample_data_files_to_results_subdirs()
        self.samples = list(self._available_samples())
        self.vcf_munger = VcfMunger()

        if len(self.samples) == 0:
            msg = 'I found no sample files (.fastq) in {}'
            raise EmptyCohortException(msg.format(self.data_dir))

        self.__vcf_stats = None

    def __repr__(self):
        tmpl = '<{} with {}>'
        return tmpl.format(self.__class__.__name__,
                           plural('sample', len(self.samples)))

    def __str__(self):
        return re.search(r'<(.*)>', self.__repr__()).groups()[0]

    def vcf_stats(self, vcf_FORMAT_field=None):
        """Create a pandas DataFrame with the values from one FORMAT field of
        the passed VCF file: 'DP' for depth, 'GQ' for genotype quality.
        The df is memoized after the first call."""
        if self.__vcf_stats is None:
            _, samples_df = VcfMunger.vcf_to_frames(self.filtered_vcf)
            self.__vcf_stats = samples_df

        if vcf_FORMAT_field:
            col_index = pd.IndexSlice[:, vcf_FORMAT_field]
            return self.__vcf_stats.loc[:, col_index]

        return self.__vcf_stats

    @staticmethod
    def read_alignment_metrics(metrics_file):
        df = pd.read_table(metrics_file, sep='\s+', comment='#')
        df['sample'] = basename(metrics_file).split('.')[0]
        # ^ Hacky, depends on the filename convention <sampleID>.extension
        return df

    # FIXME: get these metrics method into a different analyzer class
    def plot_alignment_metrics(self, metrics_files, out_path=None):
        frames = [self.read_alignment_metrics(fn) for fn in metrics_files]
        df = pd.concat(frames, ignore_index=True).dropna(axis=1)
        # Remove the info of both pairs, leave by first and second of pair.
        df_pairs = df[df['CATEGORY'] != 'PAIR']
        plotter = AlignmentMetricsPlotter(df_pairs)
        plotter.plot_and_savefig(out_path=out_path)

    def _move_sample_data_files_to_results_subdirs(self):
        glob_expr = join(self.data_dir, '*.{}'.format(Sample.reads_format))

        # Move the reads files from the data to the results dir
        # Create a folder per sample to store them.
        for reads_fp in sorted(glob(glob_expr)):
            sample_id = re.search(r'(.*).R(1|2)', basename(reads_fp)).groups(1)[0]
            sample_dir = join(self.results_dir, sample_id)
            makedirs(sample_dir, exist_ok=True)
            move(reads_fp, join(sample_dir, basename(reads_fp)))

    def _available_samples(self):
        for reads_fp in sorted(glob(join(self.results_dir, '*/*.R1.*'))):
            # Create only one Sample object per pair of reads R1-R2
            sample_dir = dirname(reads_fp)
            sample_id = basename(sample_dir)
            yield Sample(sample_id)

    def file(self, filename):
        return join(self.results_dir, filename)


class EmptyCohortException(Exception):
    pass
