import pandas as pd
import re
from glob import glob
from os.path import join, abspath, basename, expanduser
from biocodices.programs import GATK
from biocodices.helpers.language import plural
from biocodices.variant_calling import VcfMunger
from biocodices.components import Sample
from biocodices.plotters import AlignmentMetricsPlotter


class Cohort:
    def __init__(self, base_dir):
        """
        Expects the path to a directory that will have a 'data' subdirectory
        with fastq forward (R1) and reverse (R2) files from samples in it.
        """
        self.dir = abspath(expanduser(base_dir))
        self.id = basename(self.dir)
        self.data_dir = join(self.dir, 'data')
        self.results_dir = join(self.dir, 'results')
        self.samples = self._search_samples()
        self.sequencer_runs = list(set([sample.sequencer_run_id
                                        for sample in self.samples]))
        if len(self.samples) == 0:
            msg = 'I found no sample files in {}'
            raise Exception(msg.format(self.data_dir))
        self.joint_vcf = join(self.results_dir, 'joint_genotyping.vcf')

    def __repr__(self):
        tmpl = '<{}({})>'
        return tmpl.format(self.__class__.__name__, self.dir)

    def __str__(self):
        tmpl = '{} with {} from {}'
        return tmpl.format(self.__class__.__name__,
                           plural('sample', len(self.samples)),
                           ', '.join(self.sequencer_runs))

    def call_variants(self, trim_reads=True, align_reads=True,
                      create_vcfs=True, joint_genotyping=True,
                      hard_filtering=True):

        if trim_reads or align_reads or create_vcfs:
            for sample in self.samples:
                sample.call_variants(trim_reads=trim_reads,
                                     align_reads=align_reads,
                                     create_vcfs=create_vcfs)
            self.plot_alignment_metrics()

        if joint_genotyping:
            self.joint_genotyping()
            self.split_joint_vcf()

        if hard_filtering:
            for sample in self.samples:
                sample.apply_filters_to_vcf()

    def joint_genotyping(self):
        gatk = GATK()
        gvcf_list = [sample.gvcf for sample in self.samples]
        output_dir = self.results_dir
        print('\nJoint Genotyping for {}\n'.format(plural('sample',
                                                          len(self.samples))))
        self.joint_vcf = gatk.joint_genotyping(gvcf_list, output_dir)

    def split_joint_vcf(self):
        for sample in self.samples:
            outfile = sample.joint_vcf
            VcfMunger.subset(vcf=self.joint_vcf, sample_ids=[sample.id],
                             outfile=outfile)

    def plot_alignment_metrics(self):
        frames = [sample.alignment_metrics() for sample in self.samples]
        df = pd.concat(frames, ignore_index=True).dropna(axis=1)
        # Remove the info of both pairs, leave by first and second of pair.
        df_pairs = df[df['CATEGORY'] != 'PAIR']
        plotter = AlignmentMetricsPlotter(df_pairs)
        plotter.plot_and_savefig(self.dir)

    def plot_variant_calling_metrics(self):
        pass

    def _search_samples(self):
        glob_expr = join(self.data_dir, '*.{}'.format(Sample.reads_format))
        all_reads_filenames = sorted(glob(glob_expr))
        sample_ids = [re.search(r'(.*).R1', basename(reads_fn)).groups(1)[0]
                      for reads_fn in all_reads_filenames
                      if 'R1' in reads_fn]

        return [Sample(sample_id, self) for sample_id in sample_ids]
