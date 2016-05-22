from os.path import join, basename

from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller, GATK
from biocodices.helpers.general import params_dict_to_str, rename_tempfile


class VcfMunger:
    def __init__(self, sample, results_dir):
        self.sample = sample
        self.results_dir = results_dir
        self.executables = Config('executables')
        self.params = Config('parameters')
        self.gatk = GATK()

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__, self.sample.id)

    def create_snp_and_indel_vcfs(self, vcf_file):
        snps_vcf = self.gatk.select_variants(vcf_file, 'snps')
        indels_vcf = self.gatk.select_variants(vcf_file, 'indels')
        return snps_vcf, indels_vcf

    def apply_filters(self, vcf, variant_type):
        filtered_vcf = self.gatk.filter_variants_vcf(vcf, variant_type)
        return filtered_vcf
