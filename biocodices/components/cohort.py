from datetime import datetime
import re
from glob import glob
from os.path import join, abspath, basename, expanduser
import pandas as pd
from termcolor import colored

from biocodices.components import Sample
from biocodices.variant_calling import VcfMunger
from biocodices.programs import GATK
from biocodices.helpers.language import plural
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
        self.vcf_munger = VcfMunger()

        if len(self.samples) == 0:
            msg = 'I found no sample files (.fastq) in {}'
            raise Exception(msg.format(self.data_dir))

        self.unfiltered_vcf = join(self.results_dir,
                                   GATK.joint_genotyping_outfile)
        self.filtered_vcf = join(self.results_dir, GATK.hard_filtering_outfile)

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

        self.printlog('Plot some alignment metrics for the cohort.')
        self.plot_alignment_metrics()
        self.printlog('Plot the median coverage of the cohort.')
        self.plot_median_coverage()

        if joint_genotyping:
            self.printlog('Joint genotyping.')
            self.joint_genotyping()

        if hard_filtering:
            self.apply_filters_to_vcf(self.unfiltered_vcf)
            # for sample in self.samples:
                #  sample.apply_filters_to_vcf()

            self.printlog('Split the multisample VCF per sample')
            for sample in self.samples:
                self.vcf_munger.filter_samples(self.filtered_vcf, [sample.id],
                                               sample.filtered_vcf)

    def joint_genotyping(self):
        gatk = GATK()
        gvcf_list = [sample.gvcf for sample in self.samples]
        output_dir = self.results_dir
        gatk.joint_genotyping(gvcf_list, output_dir)

    def apply_filters_to_vcf(self, vcf_path):
        self.printlog('Separate SNPs and INDELs before filtering')
        variant_files = self.vcf_munger.create_snp_and_indel_vcfs(vcf_path)
        self.snps_vcf, self.indels_vcf = variant_files

        self.printlog('Apply SNP filters')
        self.snps_vcf = self.vcf_munger.apply_filters(self.snps_vcf, 'snps')

        self.printlog('Apply indel filters')
        self.indels_vcf = self.vcf_munger.apply_filters(self.indels_vcf, 'indels')

        self.printlog('Merge the filtered vcfs')
        self.vcf_munger.merge_variant_vcfs([self.snps_vcf, self.indels_vcf],
                                            outfile=self.filtered_vcf)

        #  self.printlog('Generate variant calling metrics')
        #  self.vcf_munger.generate_variant_calling_metrics(self.filtered_vcf)

    def plot_alignment_metrics(self):
        frames = [sample.read_alignment_metrics() for sample in self.samples]
        df = pd.concat(frames, ignore_index=True).dropna(axis=1)
        # Remove the info of both pairs, leave by first and second of pair.
        df_pairs = df[df['CATEGORY'] != 'PAIR']
        plotter = AlignmentMetricsPlotter(df_pairs)
        plotter.plot_and_savefig(self.results_dir)

    def median_coverages(self):
        median_coverages = {}
        for sample in self.samples:
            median_coverages[sample.id] = sample.get_median_coverage()
        median_coverages = pd.Series(median_coverages)
        return median_coverages

    #  def plot_variant_calling_metrics(self):
        #  frames = [sample.read_variant_calling_metrics()
                  #  for sample in self.samples]
        #  df = pd.concat(frames, ignore_index=True).dropna(axis=1)
        #  plotter = VariantCallingMetricsPlotter(df)
        #  plotter.plot_and_savefig(self.dir)

    def _search_samples(self):
        glob_expr = join(self.data_dir, '*.{}'.format(Sample.reads_format))
        all_reads_filenames = sorted(glob(glob_expr))
        sample_ids = [re.search(r'(.*).R1', basename(reads_fn)).groups(1)[0]
                      for reads_fn in all_reads_filenames
                      if 'R1' in reads_fn]

        return [Sample(sample_id, self) for sample_id in sample_ids]

    def printlog(self, msg):
        timestamp = datetime.now().strftime('%H:%M:%S')
        prefix = colored('[{}][{}]'.format(timestamp, self.id), 'blue')
        print('{} {}'.format(prefix, msg))
