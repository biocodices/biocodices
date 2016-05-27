from os.path import join

from biocodices.programs import AbstractGenomicsProgram, ProgramCaller


class FastQC(AbstractGenomicsProgram):
    def __init__(self):
        super(self.__class__, self).__init__('fastqc')

    def analyze_reads(self, reads_filepath, results_dir):
        command = '{} {}'.format(self.executable, reads_filepath)
        command += ' '.join(['-{} {}'.format(k, v)
                             for k, v in self.params.items()])
        command = command.format(**{'output_dir': results_dir})
        log_filepath = join(results_dir, 'fastqc')
        ProgramCaller(command).run(log_filepath=log_filepath)
