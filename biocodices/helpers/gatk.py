from os.path import dirname, join
from biocodices.helpers.config import Config
from biocodices.helpers.resource import Resource
from biocodices.helpers.program_caller import ProgramCaller


class GATK:
    def __init__(self, bam_filepath):
        self.executable = Config('executables')['GATK']
        self.reference = Resource('reference_genome')
        self.params = Config('parameters')['GATK']

        self.bam = bam_filepath
        self.realigned_bam = self.bam.replace('.bam', '.realigned.bam')
        self.recalibrated_bam = self.realigned_bam.replace('.bam',
                                                           '.recalibrated.bam')


    def realign_indels(self):
        targets_filepath = self._realigner_target_creator()
        self._indel_realigner(self.bam, targets_filepath)
        return self.realigned_bam

    def recalibrate_quality_scores(self):
        recalibration_table = self._create_recalibration_table()
        self._recalibrate_bam(self.bam, recalibration_table)
        return self.recalibrated_bam

    def _create_recalibration_table(self):
        pass

    def _realigner_target_creator(self):
        params = ['-{} {}'.format(k, v) for k, v in
                  self.params['RealignerTargetCreator'].items()]

        indels_files = [Resource('1000G_indels'), Resource('mills_indels')]
        params += ['-known {}'.format(fn) for fn in indels_files]

        targets_filepath = self.bam.replace('.bam', '.intervals')
        params_str = ' '.join(params).format(**{
            'reference_genome': self.reference,
            'input': self.bam,
            'output': targets_filepath,
            'panel_amplicons': Resource('panel_amplicons:ENPv1'),
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
            'reference_genome': self.reference,
            'input': self.bam,
            'output': self.realigned_bam,
            'target_intervals': targets_filepath,
        })

        command = '{} {}'.format(self.executable, params_str)
        log_filepath = join(dirname(self.bam), 'IndelRealigner.log')
        ProgramCaller(command).run(log_filepath=log_filepath)
