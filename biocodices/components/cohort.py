import re
from glob import glob
from os.path import join, abspath, basename, expanduser
from biocodices.programs import GATK
from biocodices.helpers.language import plural
from biocodices.variant_calling import VcfMunger
from biocodices.components import Sample


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
        self.sequencer_runs = set([sample.sequencer_run_id
                                   for sample in self.samples])
        if len(self.samples) == 0:
            msg = 'I found no sample files in {}'
            raise Exception(msg.format(self.data_dir))
        if len(self.sequencer_runs) == 1:
            self.joint_vcf = join(self.sequencer_runs[0].results_dir,
                                  'joint_genotyping.vcf')

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

        for sample in self.samples:
            sample.call_variants(trim_reads=trim_reads,
                                 align_reads=align_reads,
                                 create_vcfs=create_vcfs)

        if joint_genotyping:
            self.joint_vcf = self.joint_genotyping()

        if hard_filtering:
            for sample in self.samples:
                sample.apply_filters_to_vcf()

    def joint_genotyping(self):
        gatk = GATK()
        gvcf_list = [sample.gvcf for sample in self.samples]
        # Since the samples might come from multiple sequencer runs, we
        # pick on of them and put the joint vcf in its results directory.
        # If they're samples from only one sequence run, the path will be
        # the expected one (the results directory parent to all sample
        # directories).
        output_dir = self.sequencer_runs[0].results_dir
        print('\nJoint Genotyping for {}\n'.format(plural('sample',
                                                          len(self.samples))))
        self.joint_vcf = gatk.joint_genotyping(gvcf_list, output_dir)
        self.split_joint_vcf()

    def split_joint_vcf(self):
        for sample in self.samples:
            outfile = sample.joint_vcf
            VcfMunger.subset(vcf=self.joint_vcf, sample_ids=[sample.id],
                             outfile=outfile)

    def _search_samples(self):
        glob_expr = join(self.data_dir, '*.{}'.format(Sample.reads_format))
        all_reads_filenames = sorted(glob(glob_expr))
        sample_ids = [re.search(r'(.*).R1', basename(reads_fn)).groups(1)[0]
                      for reads_fn in all_reads_filenames
                      if 'R1' in reads_fn]

        return [Sample(sample_id, self) for sample_id in sample_ids]

