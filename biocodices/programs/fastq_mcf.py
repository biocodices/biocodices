from os.path import join, dirname

from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller
from biocodices.helpers.general import params_dict_to_str


def trim_adapters(reads_filepaths):
    """
    Expects a list of two files: forward and reverse reads of the same
    sample. It will search for an adapters file defined in resources.yml.
    """
    trimmed_fastqs = [fp.replace('.fastq', '.trimmed.fastq')
                      for fp in reads_filepaths]

    params_str = params_dict_to_str(Config.params['fastq-mcf'])
    for trimmed_fastq in trimmed_fastqs:
        params_str += ' -o {}'.format(trimmed_fastq)

    command = '{} {}'.format(Config.executables['fastq-mcf'], params_str)
    adapters_file = Resource('illumina_adapters_file')
    command += ' {} {} {}'.format(adapters_file, *reads_filepaths)

    log_filepath = join(dirname(trimmed_fastqs[0]), 'fastq-mcf')
    ProgramCaller(command).run(log_filepath=log_filepath)

    return trimmed_fastqs
