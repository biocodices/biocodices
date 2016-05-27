from os.path import dirname, join, basename
from biocodices.programs import AbstractGenomicsProgram, ProgramCaller
from biocodices.helpers.general import params_dict_to_str, rename_tempfile


class BWA(AbstractGenomicsProgram):
    def __init__(self):
        super(self.__class__, self).__init__('bwa')

    def align_to_reference(self, trimmed_reads_filepaths, ref='GRCh37'):
        """
        Expects a list of two files: forward and reverse reads of the same
        sample. It will search for a reference genome defined by Config.
        """
        outfile = trimmed_reads_filepaths[0].replace('R1.trimmed.fastq',
                                                     '.sam')
        params_str = params_dict_to_str(self.params).format(**{
            'reference_genome': self.reference_genome,
        })
        for reads_file in trimmed_reads_filepaths:
            params_str += ' {}'.format(reads_file)

        command = '{} {}'.format(self.executable, params_str)
        # redirect stdout to samfile and stderr to logfile
        log_filepath = outfile.replace('.sam', '')
        ProgramCaller(command).run(stdout_sink=outfile + '.temp',
                                   log_filepath=log_filepath)
        rename_tempfile(outfile)
