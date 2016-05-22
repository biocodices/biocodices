from os.path import join, dirname

from biocodices.helpers import Config
from biocodices.programs import ProgramCaller, GATK


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

    def merge_variant_vcfs(self, variant_vcfs, outfile):
        combined_vcf = self.gatk.combine_variants(variant_vcfs,
                                                  outfile=self.sample.vcf)
        return combined_vcf

    @staticmethod
    def subset(vcf, sample_ids, outfile):
        """
        vcf should be an absolute path to a vcf file. columns should be a
        list of the columns (i.e. samples ids) you want to subset.
        outfile should be an absolute path for the new vcf file.
        """
        command = Config('executables')['vcf-subset']
        for sample_id in sample_ids:
            command += ' -c {}'.format(sample_id)
        command += ' {}'.format(vcf)

        log_filepath = join(dirname(outfile), 'vcf-subset.log')
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
