from os.path import dirname, join
from biocodices.helpers.config import Config
from biocodices.helpers.resource import Resource
from biocodices.helpers.program_caller import ProgramCaller


class GATK:
    def __init__(self, bam_filepath):
        self.executable = Config('executables')['GATK']
        self.reference_genome = Resource('reference_genome')
        self.params = Config('parameters')['GATK']

        self.bam = bam_filepath
        self.realigned_bam = self.bam.replace('.bam', '.realigned.bam')
        self.recalibrated_bam = self.realigned_bam.replace('.bam',
                                                           '.recalibrated.bam')

    def realign_indels(self):
        targets_filepath = self._realigner_target_creator()
        self._indel_realigner(targets_filepath)
        return self.realigned_bam

    def recalibrate_quality_scores(self):
        recalibration = self._create_recalibration_table()
        self._recalibrate_bam(recalibration)
        return self.recalibrated_bam

    def _create_recalibration_table(self):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['BaseRecalibrator'].items()]

        recalibration = self.bam.replace('.bam', '.recalibration')
        indels_files = [Resource('1000G_indels'), Resource('mills_indels')]
        params += ['-knownSites {}'.format(fn) for fn in indels_files]
        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference_genome,
            'input': self.realigned_bam,
            'output': recalibration,
            'limits': Resource('panel_amplicons:ENPv1'),
            # ^ TODO: consider other panels instead of hardcoding ENPv1
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'BaseRecalibrator.log')
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
        log_filepath = join(dirname(self.bam), 'PrintReads.log')
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
            'limits': Resource('panel_amplicons:ENPv1'),
            # ^ TODO: consider other panels instead of hardcoding ENPv1
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'RealignerTargetCreator.log')
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
        log_filepath = join(dirname(self.bam), 'IndelRealigner.log')
        ProgramCaller(command).run(log_filepath=log_filepath)
