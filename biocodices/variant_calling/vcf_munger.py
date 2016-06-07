import vcf
import pandas as pd
from os.path import dirname, join

from biocodices.helpers import Config
from biocodices.programs import ProgramCaller, GATK, Picard, BcfTools


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

    def subset_samples(self, vcf_path, sample_ids, outfile):
        self.bcftools.subset_samples(vcf_path, sample_ids, outfile)
        return outfile

    def limit_regions(self, vcf_path):
        indexed_gzipped_vcf = self.bcftools.compress_and_index_vcf(vcf_path)
        return self.bcftools.limit_regions(indexed_gzipped_vcf)

    def apply_genotype_filters(self, vcf_path):
        return self.gatk.filter_genotypes(vcf_path)

    def annotate_with_snpeff(self, vcf_path):
        app_name = 'SnpEff'
        outfile = vcf_path.replace('.vcf', '.{}.vcf'.format(app_name))
        command = '{executable} {reference_genome} {in}'.format(**{
            'executable': self.executables[app_name],
            'reference_genome': self.params[app_name]['reference_genome'],
            'in': vcf_path,
        })
        log_filepath = join(dirname(vcf_path), app_name)
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
        return outfile

    def annotate_with_VEP(self, vcf_path):
        app_name = 'VEP'
        outfile = vcf_path.replace('.vcf', '.{}.vcf'.format(app_name))
        options = ['--{} {}'.format(k, v)
                   for k, v in self.params[app_name].items()]
        options_str = ' '.join(options)
        command = '{executable} {options} -i {in} -o {out}'.format(**{
            'executable': self.executables[app_name],
            'options': options_str,
            'in': vcf_path,
            'out': outfile,
        })
        log_filepath = join(dirname(vcf_path), app_name)
        ProgramCaller(command).run(log_filepath=log_filepath)
        return outfile

    @staticmethod
    def read_depth_stats_vcf(vcf_path):
        with open(vcf_path, 'r') as vcf_file:
            interval_depths = [row.INFO['IDP'] for row in vcf.Reader(vcf_file)]
        return interval_depths

    @classmethod
    def vcf_to_frames(cls, vcf_path):
        """
        Given a VCF filepath, it will return two pandas DataFrames:
            - info_df: info per variant.
            - samples_df: genotypes and genotyping data per variant / sample.
        """
        header = cls._header_from_vcf(vcf_path)
        vcf_df = pd.read_table(vcf_path, sep='\s+', comment='#', names=header,
                               index_col=['#CHROM', 'POS', 'ID'])
        info_df = vcf_df.iloc[:, 1:5]
        samples_df = cls._parse_samples_part_of_vcf(vcf_df)
        return info_df, samples_df

    @staticmethod
    def _header_from_vcf(vcf_path):
        """Auxiliary method for vcf_to_frames()"""
        with open(vcf_path, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    return line.strip().split('\t')

    @staticmethod
    def _parse_samples_part_of_vcf(vcf_df):
        """Auxiliary method for vcf_to_frames()"""
        raw_samples_df = vcf_df.iloc[:, 5:]
        assert raw_samples_df.columns[0] == 'FORMAT'
        format_col = raw_samples_df['FORMAT']

        frames = []
        for sample in raw_samples_df.columns[1:]:
            indices = list(format_col.index)[::-1]

            def genotype_dict(cell):
                ix = indices.pop()
                return dict(zip(format_col.loc[ix].split(':'),
                                cell.split(':')))

            col_as_dict = raw_samples_df[sample].map(genotype_dict).to_dict()
            df = pd.DataFrame(col_as_dict).reset_index()
            df['sample'] = sample
            df.set_index(['sample', 'index'], inplace=True)
            df = df.transpose()
            frames.append(df)

        samples_df = pd.concat(frames, axis=1).fillna('')
        return samples_df
