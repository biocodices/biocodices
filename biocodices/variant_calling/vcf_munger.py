import vcf
import pandas as pd
from os.path import dirname, join

from biocodices.helpers import Config, DB
from biocodices.programs import ProgramCaller, GATK, Picard, BcfTools


class VcfMunger:
    def __init__(self):
        self.gatk = GATK()
        self.picard = Picard()
        self.bcftools = BcfTools()

    def hard_filtering(self, vcf_path):
        """Do a hard filtering on the passed VCF. It will split it in subsets
        of INDELS and SNPS, apply different filters to each sub-VCF, and then
        merge the filtered VCFs in a new one. It returns the filepath of the
        last VCF generated."""
        snps_vcf = self.gatk.select_variants(vcf_path, 'snps')
        indels_vcf = self.gatk.select_variants(vcf_path, 'indels')
        filtered_snps_vcf = self.gatk.filter_variants_vcf(snps_vcf, 'snps')
        filtered_indels_vcf = self.gatk.filter_variants_vcf(indels_vcf,
                                                            'indels')
        merged_vcf_filepath = join(dirname(vcf_path),
                                   Config.filenames['cohort_filtered_vcf'])
        merged_vcf = self.gatk.combine_variant_vcfs(
            [filtered_snps_vcf, filtered_indels_vcf],
            outfile=merged_vcf_filepath)

        return merged_vcf

    def variant_calling_metrics(self, vcf_path):
        return self.picard.variant_calling_metrics(vcf_path)

    def limit_regions(self, vcf_path):
        """Limit the variants in a VCF within the regions defined in a .bed
        file set in the parameters config file. Generates a new VCF."""
        indexed_gzipped_vcf = self.bcftools.compress_and_index_vcf(vcf_path)
        return self.bcftools.limit_regions(indexed_gzipped_vcf)

    def annotate_with_snpeff(self, vcf_path):
        """Add SnpEff annotations in the ANN INFO field. Generates a VCF."""
        app_name = 'SnpEff'
        outfile = vcf_path.replace('.vcf', '.{}.vcf'.format(app_name))
        command = '{executable} -v {reference_genome} {in}'.format(**{
            'executable': Config.executables[app_name],
            'reference_genome': Config.params[app_name]['reference_genome'],
            'in': vcf_path,
        })
        log_filepath = join(dirname(vcf_path), app_name)
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
        return outfile

    def annotate_with_VEP(self, vcf_path):
        """Add Variant Effect Predictor's annotations in the CSQ INFO field.
        Generates a new VCF."""
        app_name = 'VEP'
        outfile = vcf_path.replace('.vcf', '.{}.vcf'.format(app_name))
        options = ['--{} {}'.format(k, v)
                   for k, v in Config.params[app_name].items()]
        options_str = ' '.join(options)
        command = '{executable} {options} -i {in} -o {out} --stats'.format(**{
            'executable': Config.executables[app_name],
            'options': options_str,
            'in': vcf_path,
            'out': outfile,
            'stats_file': vcf_path.replace('.vcf', '.VEPstats'),
        })
        log_filepath = join(dirname(vcf_path), app_name)
        ProgramCaller(command).run(log_filepath=log_filepath)
        return outfile

    def annotate_pubmed_with_DB(self, vcf_path, database, citations_table,
                                rs_column, pubmed_column):
        """Pass a database name and the name of the table with the PUBMED IDs
        per SNP. A new VCF will be written with the PUBMED IDs in a new INFO
        field."""
        citations = DB(database=database).table(citations_table)
        cited_variations = citations[rs_column]
        outfile = vcf_path.replace('.vcf', '.db.vcf')

        with open(vcf_path, 'r') as in_vcf, open(outfile, 'w') as out_vcf:
            reader = vcf.Reader(in_vcf)
            writer = vcf.Writer(out_vcf, template=reader)

            for record in reader:
                if record.ID and record.ID in cited_variations.values:
                    these_citations = citations[cited_variations == record.ID]
                    pubmed_IDs = these_citations[pubmed_column].values
                    record.INFO['PUBMED_IDs'] = ','.join(pubmed_IDs)

                writer.write_record(record)

            writer.close()

        return outfile

    # WIP:
    #  def annotate_with_cellbase(self, vcf_path):
        #  """Creates a new VCF after the passed one, adding cellbase and COSMIC
        #  data per variant in the INFO field."""
        #  outfile = vcf_path.replace('.vcf', '.db.vcf')

        #  with open(vcf_path, 'r') as in_vcf, open(outfile, 'w') as out_vcf:
            #  reader = vcf.Reader(in_vcf)
            #  writer = vcf.Writer(out_vcf, template=reader)

            #  for record in reader:
                #  # Annotate both the reference allele and the alt allele
                #  alleles = {'REF': record.REF, 'ALT': record.ALT}
                #  for which, allele in alleles.items():
                    #  kwargs = {'chromosome': record.CHROM,
                              #  'position': record.POS,
                              #  'allele': allele}
                    #  allele_info = VariantAnnotator.query_cellbase(**kwargs)
                    #  # record.INFO[which] = allele_info

    @staticmethod
    def read_depth_stats_vcf(vcf_path):
        with open(vcf_path, 'r') as vcf_file:
            interval_depths = [row.INFO['IDP'] for row in vcf.Reader(vcf_file)]
        return interval_depths

    @classmethod
    def vcf_to_frames(cls, vcf_path):
        """ Given a VCF filepath, it will return two pandas DataFrames:
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
                # FIXME: This is a nasty hack, has to be rethought
                format_at_this_row = format_col.loc[ix]
                if type(format_at_this_row) == pd.Series:
                    format_cell = format_at_this_row.iloc[0].split(':')
                else:
                    format_cell = format_at_this_row.split(':')
                sample_genotype_cell = cell.split(':')
                return dict(zip(format_cell, sample_genotype_cell))

            col_as_dict = raw_samples_df[sample].map(genotype_dict).to_dict()
            df = pd.DataFrame(col_as_dict).reset_index()
            df['sample'] = sample
            df.set_index(['sample', 'index'], inplace=True)
            df = df.transpose()
            frames.append(df)

        samples_df = pd.concat(frames, axis=1).fillna('')
        return samples_df
