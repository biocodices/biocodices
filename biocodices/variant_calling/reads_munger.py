from os.path import join, basename

from biocodices.helpers.config import Config
from biocodices.helpers.program_caller import ProgramCaller


class ReadsMunger:
    def __init__(self, results_dir):
        self.results_dir = results_dir

    def analyze_reads(self, reads_filepath, log_filepath):
        fastq_executable = Config('executables')['fastqc']
        command = '{} {} -o {}'.format(fastq_executable, reads_filepath,
                                       self.results_dir)
        ProgramCaller(command, log_filepath).run()

    def trim_adapters(self, reads_filepaths):
        trimmed_reads_filepaths = []
        for reads_filepath in reads_filepaths:
            new_fn = basename(reads_filepath).replace('.fastq', '.trimmed.fastq')
            new_filepath = join(self.results_dir, new_fn)
            trimmed_reads_filepaths.append(new_filepath)

        fastq_mcf_executable = Config('executables')['fastq-mcf']
        command = '{} -o {} -o {}'.format(fastq_mcf_executable,
                                          *trimmed_reads_filepaths)
        for k, v in Config('parameters')['fastq-mcf'].items():
            # The hacky replace is for flags with no value, i.e. -u
            command += ' -{}{}'.format(k, v)
        adapters_file = Config('resources')['illumina_adapters_file']
        command += ' {} {} {}'.format(adapters_file, *self.files['reads'])

        log_filepath = self._filepath('fastq-mcf.log')
        ProgramCaller(command, log_filepath).run()

        return trimmed_reads_filepaths
