from os.path import dirname, join
from biocodices.helpers.config import Config
from biocodices.helpers.resource import Resource
from biocodices.helpers.program_caller import ProgramCaller


class GATK:
    @classmethod
    def realign_indels(cls, bam_filepath):
        targets_filepath = cls._realigner_target_creator(bam_filepath)
        cls._indel_realigner(bam_filepath, targets_filepath)

    @staticmethod
    def _realigner_target_creator(bam_filepath):
        executable = Config('executables')['GATK'] + ' -T RealignerTargetCreator'

        params = []
        for k, v in Config('parameters')['RealignerTargetCreator'].items():
            params.append('-{} {}'.format(k, v))
        for indels_file in [Resource('1000G_indels'), Resource('mills_indels')]:
            params.append('-known {}'.format(indels_file))
        params_str = ' '.join(params).format(**{
            'reference_genome': Resource('reference_genome'),
            'input': bam_filepath,
            'panel_amplicons': Resource('panel_amplicons:ENPv1'),
            # ^ TODO: consider other panels instead of hardcoding ENPv1
            'output': bam_filepath.replace('.bam', '.realigner_targets'),
        })

        command = '{} {}'.format(executable, params_str)
        log_filepath = join(dirname(bam_filepath), 'RealignerTargetCreator.log')
        ProgramCaller(command).run(log_filepath=log_filepath)

    @staticmethod
    def _indel_realigner(bam_filepath, targets_filepath):
        pass
