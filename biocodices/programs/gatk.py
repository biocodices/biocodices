from os.path import dirname, join, basename

from biocodices.programs import AbstractGenomicsProgram, ProgramCaller
from biocodices.helpers.general import rename_tempfile


class GATK(AbstractGenomicsProgram):
    joint_genotyping_outfile = 'raw_variants.vcf'
    hard_filtering_outfile = 'filtered.vcf'

    def __init__(self):
        super(self.__class__, self).__init__('GATK')

    def run(self, module_name, infile, outfile, extra_params_variables={},
            extra_output_extension=None, log_label=None, task_subtype=None,
            extra_params_str=None):
        params_dict = self.params[module_name]
        if task_subtype:
            params_dict = params_dict[task_subtype]
            log_label = '{}_{}'.format(module_name, task_subtype)
        params = ['-{} {}'.format(k, v) for k, v in params_dict.items()]
        params_variables = {
            'module_name': module_name,
            'reference_genome': self.reference_genome,
            'known_variants': self.known_variants,
            'limits': self.panel_amplicons,
            'input': infile,
            'output': outfile + '.temp',
        }
        params_variables.update(extra_params_variables)
        params_str = ' '.join(params).format(**params_variables)
        if extra_params_str:
            params_str += (' ' + extra_params_str)
        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(infile or outfile),
                            (log_label or module_name))
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, extra_output_extension)

    def realign_reads_around_indels(self, bam):
        targets_filepath = self._realigner_target_creator(bam)
        realigned_bam = self._indel_realigner(targets_filepath, bam)
        return realigned_bam

    def recalibrate_quality_scores(self, realigned_bam):
        recalibration = self._create_recalibration_table(realigned_bam)
        recalibrated_bam = self._recalibrate_bam(realigned_bam, recalibration)
        return recalibrated_bam

    def create_depth_vcf(self, recalibrated_bam):
        """Expects a recalibrated BAM and outputs an depth stats VCF"""
        outfile = recalibrated_bam.replace('.bam', '.depth_stats.vcf')
        self.run('DiagnoseTargets', recalibrated_bam, outfile)
        return outfile

    def create_vcf(self, recalibrated_bam):
        """
        Expects a recalibrated BAM as infile, ready for the variant calling.
        """
        outfile = recalibrated_bam.replace('.bam', '.vcf')
        self.run('HaplotypeCaller', recalibrated_bam, outfile,
                 extra_output_extension='idx', task_subtype='vcf')
        return outfile

    def create_gvcf(self, recalibrated_bam):
        """
        Expects a recalibrated BAM as infile, ready for the variant calling.
        """
        outfile = recalibrated_bam.replace('.bam', '.g.vcf')
        self.run('HaplotypeCaller', recalibrated_bam, outfile,
                 extra_output_extension='idx', task_subtype='gvcf')
        return outfile

    def joint_genotyping(self, gvcf_list, output_dir):
        outfile = join(output_dir, self.joint_genotyping_outfile)
        params_str = ''
        for gvcf_filename in gvcf_list:
            params_str += ' --variant {}'.format(gvcf_filename)

        # The None is bc joint genotyping doesn't take an I (input) parameter
        self.run('GenotypeGVCFs', None, outfile, extra_params_str=params_str,
                 extra_output_extension='idx')

        return outfile

    def select_variants(self, vcf, variant_type):
        variant_vcf = vcf.replace('.vcf', '.{}.vcf'.format(variant_type))
        module_name = 'SelectVariants'
        log_label = '{}_{}'.format(module_name, variant_type)
        self.run(module_name, vcf, variant_vcf, task_subtype=variant_type,
                 extra_output_extension='idx', log_label=log_label)
        return variant_vcf

    def filter_variants_vcf(self, vcf, variant_type):
        """
        Apply different filters listed in the parameters.yml config file to
        a vcf splitted earlier by SNP vs INDEL.
        """
        module_name = 'VariantFiltration'
        # The filters should be specific for the kind of variant being
        # filtered. Applying them in order is important for a consisting
        # naming of the output files throughout different runs.

        input_vcf = vcf
        outfile = input_vcf.replace('.vcf', '.filtered.vcf')
        filters = self.params['{}_filters'.format(module_name)][variant_type]
        filter_name = '{}_filter'.format(variant_type)
        filter_expression = ' || '.join(filters)
        params_str = '--filterName {} --filterExpression "{}"'.format(
            filter_name, filter_expression)
        log_label = '{}_{}'.format(module_name, filter_name)

        self.run(module_name, vcf, outfile, log_label=log_label,
                 extra_params_str=params_str, extra_output_extension='idx')
        return outfile

    def filter_genotypes(self, vcf):
        module_name = 'VariantFiltration'
        filter_name = 'genotype'
        outfile = vcf.replace('.vcf', '.geno-filtered.vcf')
        filters = self.params['{}_filters'.format(module_name)][filter_name]
        filter_expression = ' || '.join(filters)
        params_str = '--genotypeFilterName {} '.format(filter_name)
        params_str += '--genotypeFilterExpression "{}"'.format(filter_expression)
        log_label = '{}_{}'.format(module_name, filter_name)

        self.run(module_name, vcf, outfile, log_label=log_label,
                 extra_params_str=params_str, extra_output_extension='idx')
        return outfile

    def combine_variant_vcfs(self, variant_vcfs, outfile):
        """
        Written to merge INDEL and SNP vcfs from the same sample.
        """
        # Check the VCFs are from the same sample (just checks the filename).
        sample_ids = [basename(vcf).split('.')[0] for vcf in variant_vcfs]
        sample_ids = list(set(sample_ids))
        if len(sample_ids) > 1:
            msg = 'Are you trying to merge vcfs from different samples? {}'
            raise Exception(msg.format(variant_vcfs))

        params_str = ''
        # GATK CombineVariants needs to tag the files being merged to then
        # assign them a priority in the merging process. I just use a dummy
        # tag since every variant will be included (I guess there will be no
        # conflicts in the merged data?).
        tag = sample_ids[0]
        for variant_vcf in variant_vcfs:
            params_str += ' -V:{} {}'.format(tag, variant_vcf)
        dummy_tags = [tag for vcf in variant_vcfs]
        params_str += ' -priority {},{}'.format(*dummy_tags)

        # infile parameter is None because the actual inputs will be passed
        # with the -V <variant> flag, leaving the -I flag unused.
        self.run('CombineVariants', None, outfile, extra_params_str=params_str,
                 extra_output_extension='idx')
        return outfile

    def _realigner_target_creator(self, bam):
        target_intervals = bam.replace('.bam', '.intervals')
        extra_params_str = ''
        for indels_file in self.known_indels:
            extra_params_str += ' -known {}'.format(indels_file)

        self.run('RealignerTargetCreator', bam, target_intervals,
                 extra_params_str=extra_params_str)
        return target_intervals

    def _indel_realigner(self, targets_filepath, bam):
        realigned_bam = bam.replace('.bam', '.realigned.bam')
        extra_params_variables = {
            'target_intervals': targets_filepath,
        }
        self.run('IndelRealigner', bam, realigned_bam,
                 extra_output_extension='bai',
                 extra_params_variables=extra_params_variables)
        return realigned_bam

    def _create_recalibration_table(self, realigned_bam):
        recalibration = realigned_bam.replace('.bam', '.recalibration')
        extra_params_str = ''
        for indels_file in self.known_indels:
            extra_params_str += ' -knownSites {}'.format(indels_file)

        self.run('BaseRecalibrator', realigned_bam, recalibration,
                 extra_params_str=extra_params_str)
        return recalibration

    def _recalibrate_bam(self, realigned_bam, recalibration):
        recalibrated_bam = recalibration.replace('.recalibration',
                                                 '.recalibrated.bam')
        extra_params_variables = {
            'recalibration': recalibration,
        }
        self.run('PrintReads', realigned_bam, recalibrated_bam,
                 extra_params_variables=extra_params_variables,
                 extra_output_extension='bai')
        return recalibrated_bam
