from os.path import join, basename

from biocodices.helpers.config import Config
from biocodices.helpers.program_caller import ProgramCaller


class ReadsMunger:
    def __init__(self, results_dir):
        self.results_dir = results_dir
        self.fastq_executable = Config('executables')['fastqc']
        self.fastq_mcf_executable = Config('executables')['fastq-mcf']

    def analyze_reads(self, reads_filepath):
        command = '{} {} -o {}'.format(self.fastq_executable, reads_filepath,
                                       self.results_dir)
        ProgramCaller(command, self._log_filepath('fastqc')).run()

    def trim_adapters(self, reads_filepaths):
        """
        Expects two files: forward and reverse reads of the same sample.
        It will use an adapters file whose location should be defined previously.
        """
        trimmed_reads_filepaths = []
        for reads_filepath in reads_filepaths:
            new_fn = basename(reads_filepath).replace('.fastq', '.trimmed.fastq')
            new_filepath = join(self.results_dir, new_fn)
            trimmed_reads_filepaths.append(new_filepath)

        command = '{} -o {} -o {}'.format(self.fastq_mcf_executable,
                                          *trimmed_reads_filepaths)
        for k, v in Config('parameters')['fastq-mcf'].items():
            command += ' -{}{}'.format(k, v)
        adapters_file = Config('resources')['illumina_adapters_file']
        command += ' {} {} {}'.format(adapters_file, *reads_filepaths)
        print(command)

        log_filepath = self._log_filepath('fastq-mcf')
        ProgramCaller(command, log_filepath).run()

        return trimmed_reads_filepaths

    def _log_filepath(self, label):
        return join(self.results_dir, label + '.log')
