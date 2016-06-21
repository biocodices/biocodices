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

    def hard_filtering(self, vcf_path, out_path):
        """Do a hard filtering on the passed VCF. It will split it in subsets
        of INDELS and SNPS, apply different filters to each sub-VCF, and then
        merge the filtered VCFs in a new one at <out_path>."""
        raw_snps_vcf = self.gatk.select_variants(vcf_path, 'snps')
        raw_indels_vcf = self.gatk.select_variants(vcf_path, 'indels')
        snps_vcf = self.gatk.filter_variants_vcf(raw_snps_vcf, 'snps')
        indels_vcf = self.gatk.filter_variants_vcf(raw_indels_vcf, 'indels')
        return self.gatk.combine_variant_vcfs([snps_vcf, indels_vcf],
                                              outfile=out_path)

    def realign_reads_around_indels(self, bam):
        targets_filepath = self.gatk.realigner_target_creator(bam)
        realigned_bam = self.gatk.indel_realigner(targets_filepath, bam)
        return realigned_bam

    def recalibrate_quality_scores(self, realigned_bam, out_path=None):
        recalibration = self.gatk.create_recalibration_table(realigned_bam)
        recalibrated_bam = self.gatk.recalibrate_bam(realigned_bam, recalibration,
                                                     out_path=out_path)
        return recalibrated_bam

    @classmethod
    def get_median_coverage(cls, depth_vcf):
        depth_stats = cls.read_depth_stats_vcf(depth_vcf)
        return pd.Series(depth_stats).median()

    def limit_regions(self, vcf_path, out_path=None):
        """Limit the variants in a VCF within the regions defined in a .bed
        file set in the parameters config file. Generates a new VCF."""
        indexed_gzipped_vcf = self.bcftools.compress_and_index_vcf(vcf_path)
        return self.bcftools.limit_regions(indexed_gzipped_vcf, out_path)

    def annotate_with_snpeff(self, vcf_path, out_path=None):
        """Add SnpEff annotations in the ANN INFO field. Generates a VCF."""
        app_name = 'SnpEff'
        outfile = out_path or vcf_path.replace('.vcf', '.{}.vcf'.format(app_name))
        command = '{executable} -v {reference_genome} {in}'.format(**{
            'executable': Config.executables[app_name],
            'reference_genome': Config.params[app_name]['reference_genome'],
            'in': vcf_path,
        })
        log_filepath = join(dirname(vcf_path), app_name)
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
        return outfile

    def annotate_with_VEP(self, vcf_path, out_path=None):
        """Add Variant Effect Predictor's annotations in the CSQ INFO field.
        Generates a new VCF."""
        app_name = 'VEP'
        outfile = out_path or vcf_path.replace('.vcf', '.{}.vcf'.format(app_name))
        options = ['--{} {}'.format(k, v)
                   for k, v in Config.params[app_name].items()]
        options_str = ' '.join(options).format(**{
            'stats_file': vcf_path.replace('.vcf', '.VEP_stats.html'),
        })
        command = '{executable} {options} -i {in} -o {out}'.format(**{
            'executable': Config.executables[app_name],
            'options': options_str,
            'in': vcf_path,
            'out': outfile,
        })
        log_filepath = join(dirname(vcf_path), app_name)
        ProgramCaller(command).run(log_filepath=log_filepath)
        return outfile

    def annotate_pubmed_with_DB(self, vcf_path, database, citations_table,
                                rs_column, pubmed_column, out_path=None):
        """Pass a database name and the name of the table with the PUBMED IDs
        per SNP. A new VCF will be written with the PUBMED IDs in a new INFO
        field."""
        citations = DB(database=database).table(citations_table)
        cited_variations = citations[rs_column]
        outfile = out_path or vcf_path.replace('.vcf', '.db.vcf')

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
        raw_samples_df = vcf_df.iloc[:, 5:]
        samples_df = cls._parse_samples_part_of_vcf(raw_samples_df)
        return info_df, samples_df

    @staticmethod
    def _header_from_vcf(vcf_path):
        """Auxiliary method for vcf_to_frames()"""
        with open(vcf_path, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    return line.strip().split('\t')

    @staticmethod
    def _parse_samples_part_of_vcf(raw_samples_df):
        """Auxiliary method for vcf_to_frames()"""

        format_series = raw_samples_df['FORMAT'].to_dict()
        # ^ Converting to_dict() removes rows with a duplicated index.

        new_df = pd.DataFrame({})

        for sample_id, gt_column in raw_samples_df.iloc[:, 1:].iteritems():
            for ix, gt_data in gt_column.str.split(':').iteritems():
                format_keys = format_series[ix].split(':')
                row_as_dict = {sample_id: dict(zip(format_keys, gt_data))}

                # unstack() creates a multi-index SAMPLE_ID > [GT, AD, DP, GQ]
                new_series = pd.DataFrame(row_as_dict).unstack()
                new_df[ix, sample_id] = new_series

        new_df = new_df.transpose()
        new_df.index = raw_samples_df.index.drop_duplicates()
        # ^ Restores the multi-index because at this point it was
        #   converted to a simple index of tuples.

        return new_df
