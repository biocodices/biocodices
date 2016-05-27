from os.path import dirname, join, basename
from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller
from biocodices.helpers.general import params_dict_to_str, rename_tempfile


class GATK:
    def __init__(self):
        self.executable = Config('executables')['GATK']
        self.params = Config('parameters')['GATK']
        self.reference_genome = Resource('reference_genome')
        self.known_indels = [Resource('indels:1000G'),
                             Resource('indels:mills')]
        self.panel_amplicons = Resource('panel_amplicons')
        self.known_variants = Resource('dbsnp:GRCh37')

    def run(self, module_name, infile, outfile, extra_params_variables={},
            extra_output_extension=None, log_label=None, task_subtype=None,
            extra_params_str=None):
        params_dict = self.params[module_name]
        if module_name == 'HaplotypeCaller':
            params_dict = params_dict[task_subtype]
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
            params_str += extra_params_str
        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(infile), (log_label or module_name))
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

    def create_vcf(self, recalibrated_bam):
        """
        Expects a recalibrated BAM as infile, ready for the variant calling.
        """
        outfile = recalibrated_bam.replace('.bam', '.vcf')
        self.run('HaplotypeCaller', recalibrated_bam, outfile,
                 extra_output_extension='idx', task_subtype='vcf',
                 log_label='HaplotypeCaller_vcf')

    def create_gvcf(self, recalibrated_bam):
        """
        Expects a recalibrated BAM as infile, ready for the variant calling.
        """
        outfile = recalibrated_bam.replace('.bam', '.g.vcf')
        self.run('HaplotypeCaller', recalibrated_bam, outfile,
                 extra_output_extension='idx', task_subtype='gvcf',
                 log_label='HaplotypeCaller_gvcf')

    def joint_genotyping(self, gvcf_list, output_dir):
        _ = None  # dummy 'infile' argument for run()
        outfile = join(output_dir, 'joint_genotyping.vcf')
        params_str = ''
        for gvcf_filename in gvcf_list:
            params_str += ' --variant {}'.format(gvcf_filename)

        self.run('GenotypeGVCFs', _, outfile, extra_params_str=params_str,
                 extra_output_extension='idx')

    def select_variants(self, vcf, variant_type):
        variant_vcf = vcf.replace('.vcf', '.{}.vcf'.format(variant_type))
        module_name = 'SelectVariants'
        log_label = '{}_{}'.format(module_name, variant_type)
        self.run(module_name, vcf, variant_vcf,
                 extra_output_extension='idx', log_label=log_label)
        return variant_vcf

    def filter_variants_vcf(self, vcf, variant_type):
        """
        Apply different filters listed in the parameters.yml config file to
        a vcf splitted earlier by SNP vs INDEL.
        """
        module_name = 'VariantFiltration'
        filters = self.params['{}_filters'.format(module_name)]
        filter_order = filters['{}_order'.format(variant_type)]
        # The filters should be specific for the kind of variant being
        # filtered. Applying them in order is important for a consisting
        # naming of the output files throughout different runs.

        input_vcf = vcf
        for filter_name in filter_order:
            outfile = input_vcf.replace('.vcf', '.{}.vcf'.format(filter_name))
            filter_expression = filters[variant_type][filter_name]
            params_str = ' --filterName {} --filterExpression {}'
            params_str = params_str.format(filter_name, filter_expression)
            log_label = '{}_{}_{}'.format(module_name, variant_type,
                                          filter_name)
            self.run(module_name, vcf, outfile, log_label=log_label,
                     extra_params_str=params_str, extra_output_extension='idx')
            input_vcf = outfile
            # The filtered vcf will be input for the next filtering round.

        return outfile  # last vcf file created, with all filters applied

    def combine_variants(self, variant_vcfs, outfile):
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
        tag = 'DummyTag'
        for variant_vcf in variant_vcfs:
            params_str += ' -V:{} {}'.format(tag, variant_vcf)
        dummy_tags = [tag for vcf in variant_vcfs]
        params_str += ' -priority {},{}'.format(*dummy_tags)

        # infile parameter is None because the actual inputs will be passed
        # with the -V <variant> flag, leaving the -I flag unused.
        self.run('CombineVariants', None, outfile, extra_params_str=params_str)
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
