from os.path import join, dirname
import vcf

from biocodices.helpers import Config
from biocodices.programs import ProgramCaller, GATK, Picard  #, VcfTools


class VcfMunger:
    def __init__(self, sample, results_dir):
        self.sample = sample
        self.results_dir = results_dir
        self.executables = Config('executables')
        self.params = Config('parameters')
        self.gatk = GATK()
        self.picard = Picard()
        # self.vcf_tools = VcfTools()

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__, self.sample.id)

    def create_snp_and_indel_vcfs(self, vcf_file):
        snps_vcf = self.gatk.select_variants(vcf_file, 'snps')
        indels_vcf = self.gatk.select_variants(vcf_file, 'indels')
        return snps_vcf, indels_vcf

    def apply_filters(self, vcf_path, variant_type):
        return self.gatk.filter_variants_vcf(vcf_path, variant_type)

    def merge_variant_vcfs(self, variant_vcfs, outfile):
        return self.gatk.combine_variant_vcfs(variant_vcfs, outfile=outfile)

    def variant_calling_metrics(self, vcf_path):
        return self.picard.variant_calling_metrics(vcf_path)

    #  def remove_variants_outside_limits(self, vcf):
        #  return self.vcf_tools.remove_variants_outside_limits(vcf)

    @staticmethod
    def read_depth_stats_vcf(vcf_path):
        with open(vcf_path, 'r') as vcf_file:
            interval_depths = [row.INFO['IDP'] for row in vcf.Reader(vcf_file)]
        return interval_depths

    @staticmethod
    def subset(vcf_path, sample_ids, outfile):
        """
        vcf should be an absolute path to a [g]VCF file. columns should be a
        list of the columns (i.e. samples ids) you want to subset.
        outfile should be an absolute path for the new vcf file.
        """
        command = Config('executables')['vcf-subset']
        for sample_id in sample_ids:
            command += ' -c {}'.format(sample_id)
        command += ' {}'.format(vcf_path)

        log_filepath = join(dirname(outfile), 'vcf-subset.log')
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
