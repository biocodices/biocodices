from os.path import dirname
from collections import OrderedDict
import vcf
import pandas as pd

from biocodices.helpers import Config
from biocodices.programs import GATK, Picard, BcfTools


class VcfMunger:
    def __init__(self):
        self.executables = Config('executables')
        self.params = Config('parameters')
        self.gatk = GATK()
        self.picard = Picard()
        self.bcftools = BcfTools()

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

    def filter_samples(self, vcf_path, sample_ids, outfile):
        self.bcftools.filter_samples(vcf_path, sample_ids, outfile)
        return outfile

    @staticmethod
    def read_depth_stats_vcf(vcf_path):
        with open(vcf_path, 'r') as vcf_file:
            interval_depths = [row.INFO['IDP'] for row in vcf.Reader(vcf_file)]
        return interval_depths

    @staticmethod
    def vcf_to_frames(vcf_path):
        """
        Given a VCF filepath, it will return two pandas DataFrames:
            - info_df: info per variant.
            - samples_df: genotypes and genotyping data per variant / sample.
        """
        info_df = pd.DataFrame({})
        samples_df = pd.DataFrame({})

        with open(vcf_path, 'r') as vcf_file:
            for record in vcf.Reader(vcf_file):
                entry = OrderedDict()
                entry['chrom'] = int(record.CHROM)
                entry['pos'] = record.POS
                entry['rs_id'] = record.ID
                entry['ref'] = record.REF
                entry['alt'] = record.ALT
                entry['qual'] = record.QUAL
                entry['filter'] = record.FILTER
                entry['format'] = record.FORMAT

                format_keys = record.FORMAT.split(':')

                for key, value in record.INFO.items():
                    entry[key] = value

                info_df = info_df.append(pd.Series(entry), ignore_index=True)

                samples_entry = OrderedDict()
                for sample in record.samples:
                    samples_entry[sample.sample] = {}
                    for key in format_keys:
                        samples_entry[sample.sample][key] = sample[key]

                df = pd.DataFrame(samples_entry).unstack().to_frame().transpose()
                df['rs_id'] = record.ID
                df['chrom'] = record.CHROM
                df['pos'] = record.POS
                samples_df = samples_df.append(df)

        info_df = info_df[info_df.columns[::-1]]
        info_df.set_index(['chrom', 'pos', 'rs_id'], inplace=True)
        samples_df = samples_df.fillna('').set_index(['chrom', 'pos', 'rs_id'])

        return info_df, samples_df
