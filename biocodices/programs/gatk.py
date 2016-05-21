from os.path import dirname, join
from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller
from biocodices.helpers.general import params_dict_to_str


class GATK:
    def __init__(self):
        self.executable = Config('executables')['GATK']
        self.reference_genome = Resource('reference_genome')
        self.params = Config('parameters')['GATK']

    def set_bamfile(self, bam_filepath):
        # TODO: this doesn't convince me. Think of a better way.
        self.bam = bam_filepath
        self.realigned_bam = self.bam.replace('.bam', '.realigned.bam')
        self.recalibrated_bam = self.realigned_bam.replace('.bam',
                                                           '.recalibrated.bam')
        self.vcf = self.bam.replace('.bam', '.vcf')
        self.gvcf = self.bam.replace('.bam', '.g.vcf')

    def realign_indels(self):
        targets_filepath = self._realigner_target_creator()
        self._indel_realigner(targets_filepath)
        return self.realigned_bam

    def recalibrate_quality_scores(self):
        recalibration = self._create_recalibration_table()
        self._recalibrate_bam(recalibration)
        return self.recalibrated_bam

    def create_vcf(self):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['HaplotypeCaller']['vcf'].items()]
        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'limits': Resource('panel_amplicons'),
            'input': self.recalibrated_bam,
            'output': self.vcf,
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'HaplotypeCaller_vcf')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def create_gvcf(self):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['HaplotypeCaller']['gvcf'].items()]
        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'limits': Resource('panel_amplicons'),
            'input': self.recalibrated_bam,
            'output': self.gvcf,
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'HaplotypeCaller_gvcf')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def joint_genotyping(self, gvcf_list, output_gvcf_filepath):
        params_dict = self.params['GATK']['GenotypeGVCFs']
        params_str = params_dict_to_str(params_dict)
        params_str = params_str.format(**{
            'reference_genome': self.reference_genome,
            'reference_variants': self.reference_variants,
            'output': output_gvcf_filepath,
        })
        for gvcf_filename in gvcf_list:
            params_str += ' --variant {}'.format(gvcf_filename)

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(output_gvcf_filepath),
                            'GenotypeGVCFs')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def _create_recalibration_table(self):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['BaseRecalibrator'].items()]

        recalibration = self.bam.replace('.bam', '.recalibration')
        indels_files = [Resource('1000G_indels'), Resource('mills_indels')]
        params += ['-knownSites {}'.format(fn) for fn in indels_files]
        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'limits': Resource('panel_amplicons'),
            'input': self.realigned_bam,
            'output': recalibration,
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'BaseRecalibrator')
        ProgramCaller(command).run(log_filepath=log_filepath)

        return recalibration

    def _recalibrate_bam(self, recalibration):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['PrintReads'].items()]

        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'input': self.realigned_bam,
            'recalibration': recalibration,
            'output': self.recalibrated_bam,
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'PrintReads')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def _realigner_target_creator(self):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['RealignerTargetCreator'].items()]

        indels_files = [Resource('1000G_indels'), Resource('mills_indels')]
        params += ['-known {}'.format(fn) for fn in indels_files]

        targets_filepath = self.bam.replace('.bam', '.intervals')

        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'input': self.bam,
            'output': targets_filepath,
            'limits': Resource('panel_amplicons'),
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'RealignerTargetCreator')
        ProgramCaller(command).run(log_filepath=log_filepath)

        return targets_filepath

    def _indel_realigner(self, targets_filepath):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['IndelRealigner'].items()]
        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'input': self.bam,
            'output': self.realigned_bam,
            'target_intervals': targets_filepath,
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'IndelRealigner')
        ProgramCaller(command).run(log_filepath=log_filepath)