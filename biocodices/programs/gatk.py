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

    def set_bamfile(self, bam_filepath):
        # TODO: this doesn't convince me. Think of a better way.
        self.bam = bam_filepath
        self.realigned_bam = self.bam.replace('.bam', '.realigned.bam')
        self.recalibrated_bam = self.realigned_bam.replace('.bam',
                                                           '.recalibrated.bam')
        self.vcf = self.bam.replace('.bam', '.vcf')
        self.gvcf = self.bam.replace('.bam', '.g.vcf')

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

    def create_vcf(self, infile):
        """
        Expects a recalibrated BAM as infile, ready for the variant calling.
        """
        outfile = infile.replace('.bam', '.vcf')
        self.run('HaplotypeCaller', infile, outfile,
                 extra_output_extension='bai', task_subtype='vcf',
                 log_label='HaplotypeCaller_vcf')

    def create_gvcf(self, infile):
        """
        Expects a recalibrated BAM as infile, ready for the variant calling.
        """
        outfile = infile.replace('.bam', '.g.vcf')
        self.run('HaplotypeCaller', infile, outfile,
                 extra_output_extension='bai', task_subtype='gvcf',
                 log_label='HaplotypeCaller_gvcf')

    def realign_indels(self):
        targets_filepath = self._realigner_target_creator()
        self._indel_realigner(targets_filepath)
        return self.realigned_bam

    def recalibrate_quality_scores(self):
        recalibration = self._create_recalibration_table()
        self._recalibrate_bam(recalibration)
        return self.recalibrated_bam

    def joint_genotyping(self, gvcf_list, output_dir):
        _ = None  # dummy 'infile' argument for run()
        outfile = join(output_dir, 'joint_genotyping.vcf')
        params_str = ''
        for gvcf_filename in gvcf_list:
            params_str += ' --variant {}'.format(gvcf_filename)

        self.run('GenotypeGVCFs', _, outfile, extra_params_str=params_str,
                 extra_output_extension='idx')

    def joint_genotyping(self, gvcf_list, output_dir):
        outfile = join(output_dir, 'joint_genotyping.vcf')
        params_dict = self.params['GenotypeGVCFs']
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'known_variants': self.known_variants,
            'output': outfile + '.temp',
        })
        for gvcf_filename in gvcf_list:
            params_str += ' --variant {}'.format(gvcf_filename)

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(output_dir, 'GenotypeGVCFs')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, 'idx')
        return outfile

    def select_variants(self, vcf, variant_type):
        outfile = vcf.replace('.vcf', '.{}.vcf'.format(variant_type))
        params_dict = self.params['SelectVariants'][variant_type]
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'input': vcf,
            'output': outfile + '.temp',
        })
        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(vcf),
                            'SelectVariants_{}'.format(variant_type))
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, 'idx')
        return outfile

    def filter_variants_vcf(self, vcf, variant_type):
        filters = self.params['VariantFiltration_filters']
        filter_order = filters['{}_order'.format(variant_type)]
        params_dict = self.params['VariantFiltration']
        input_vcf = vcf
        for filter_name in filter_order:
            filter_expression = filters[variant_type][filter_name]
            outfile = input_vcf.replace('.vcf', '.{}.vcf'.format(filter_name))
            params_str = params_dict_to_str(params_dict).format(**{
                'reference_genome': self.reference_genome,
                'input': input_vcf,
                'output': outfile + '.temp',
            })
            params_str += ' --filterName {}'.format(filter_name)
            params_str += ' --filterExpression "{}"'.format(filter_expression)
            command = '{} {}'.format(self.executable, params_str)
            log_filepath = join(dirname(input_vcf),
                                'VariantFiltration_{}_{}'.format(variant_type,
                                                                 filter_name))
            ProgramCaller(command).run(log_filepath=log_filepath)
            rename_tempfile(outfile, 'idx')
            input_vcf = outfile
            # The filtered vcf will be input for the next filtering round.

        return outfile  # last vcf file created, with all filters applied

    def combine_variants(self, variant_vcfs, outfile):
        sample_ids = [basename(vcf).split('.')[0] for vcf in variant_vcfs]
        sample_ids = list(set(sample_ids))
        if len(sample_ids) > 1:
            msg = 'Are you trying to merge vcfs from different samples? {}'
            raise Exception(msg.format(variant_vcfs))
        params_dict = self.params['CombineVariants']
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'output': outfile + '.temp',
        })

        # GATK CombineVariants needs to tag the files being merged to then
        # assign them a priority in the merging process. I just use a dummy
        # tag since every variant will be included (I guess there will be no
        # conflicts in the merged data?).
        tag = 'DummyTag'
        for variant_vcf in variant_vcfs:
            params_str += ' -V:{} {}'.format(tag, variant_vcf)
        dummy_tags = [tag for vcf in variant_vcfs]
        params_str += ' -priority {},{}'.format(*dummy_tags)

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(outfile), 'CombineVariants')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile)
        return outfile

    def _create_recalibration_table(self):
        outfile = self.bam.replace('.bam', '.recalibration')
        params_dict = self.params['BaseRecalibrator']
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'limits': Resource('panel_amplicons'),
            'input': self.realigned_bam,
            'output': outfile + '.temp',
        })
        for indels_file in self.known_indels:
            params_str += ' -knownSites {}'.format(indels_file)

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'BaseRecalibrator')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile)

        return outfile

    def _recalibrate_bam(self, recalibration):
        outfile = self.recalibrated_bam
        params_dict = self.params['PrintReads']
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'input': self.realigned_bam,
            'recalibration': recalibration,
            'output': outfile + '.temp',
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'PrintReads')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, 'bai')

    def _realigner_target_creator(self):
        targets_filepath = self.bam.replace('.bam', '.intervals')
        outfile = targets_filepath

        params_dict = self.params['RealignerTargetCreator']
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'input': self.bam,
            'limits': Resource('panel_amplicons'),
            'output': outfile + '.temp',
        })
        for indels_file in self.known_indels:
            params_str += ' -known {}'.format(indels_file)

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'RealignerTargetCreator')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile)

        return targets_filepath

    def _indel_realigner(self, targets_filepath):
        outfile = self.realigned_bam
        params_dict = self.params['IndelRealigner']
        params_str = params_dict_to_str(params_dict).format(**{
            'reference_genome': self.reference_genome,
            'input': self.bam,
            'target_intervals': targets_filepath,
            'output': outfile + '.temp',
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'IndelRealigner')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, 'bai')
