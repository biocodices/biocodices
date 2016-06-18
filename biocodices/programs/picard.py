# from os import rename
from os.path import dirname, join
from biocodices.programs import AbstractGenomicsProgram, ProgramCaller
from biocodices.helpers.general import rename_tempfile


class Picard(AbstractGenomicsProgram):
    def __init__(self):
        super(self.__class__, self).__init__('picard-tools')

    def run(self, module_name, infile, outfile, extra_params_variables={},
            extra_output_extension=None, out_rename=True):
        params_dict = self.params[module_name]
        params = ['{}={}'.format(k, v) for k, v in params_dict.items()]
        params_variables = {
            'reference_genome': self.reference_genome,
            'reference_genome_dict': self.reference_genome_dict,
            'known_variants': self.known_variants,
            'input': infile,
            'output': outfile + '.temp',
        }
        params_variables.update(extra_params_variables)
        params_str = ' '.join(params).format(**params_variables)
        command = '{} {} {}'.format(self.executable, module_name, params_str)
        log_filepath = join(dirname(infile), module_name)
        ProgramCaller(command).run(log_filepath=log_filepath)
        if out_rename:
            rename_tempfile(outfile, extra_output_extension)

    def alignment_metrics(self, recalibrated_bam):
        outfile = recalibrated_bam.replace('.bam', '.alignment_metrics.tsv')
        self.run('CollectAlignmentSummaryMetrics', recalibrated_bam, outfile)
        return outfile

    def add_or_replace_read_groups(self, sam_path, sample_id,
                                   sample_library_id, sequencer_run_id):
        outfile = sam_path.replace('.sam', '.bam')
        extra_params_variables = {
            'sample_id': sample_id,
            'library_id': sample_library_id,
            'ngs_id': sequencer_run_id,
        }
        self.run('AddOrReplaceReadGroups', sam_path, outfile,
                 extra_params_variables, 'bai')
        return outfile
