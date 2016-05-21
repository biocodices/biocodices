from os.path import join, basename

from biocodices.helpers.config import Config
from biocodices.helpers.resource import Resource
from biocodices.helpers.program_caller import ProgramCaller
from biocodices.helpers.gatk import GATK


class ReadsMunger:
    def __init__(self, sample, results_dir):
        self.sample = sample
        self.results_dir = results_dir
        self.executables = Config('executables')
        self.params = Config('parameters')

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__, self.sample.id)

    def analyze_reads(self, reads_filepath):
        command = '{} {} -o {}'.format(self.executables['fastqc'],
                                       reads_filepath, self.results_dir)
        log_filepath = self._log_filepath('fastqc')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def trim_adapters(self, reads_filepaths):
        """
        Expects a list of two files: forward and reverse reads of the same
        sample. It will search for an adapters file defined by Config.
        """
        trimmed_reads_filepaths = []
        for reads_filepath in reads_filepaths:
            new_fn = basename(reads_filepath).replace('.fastq', '.trimmed.fastq')
            new_filepath = join(self.results_dir, new_fn)
            trimmed_reads_filepaths.append(new_filepath)

        command = '{} -o {} -o {}'.format(self.executables['fastq-mcf'],
                                          *trimmed_reads_filepaths)
        for k, v in self.params['fastq-mcf'].items():
            command += ' -{}{}'.format(k, v)
        adapters_file = Resource('illumina_adapters_file')
        command += ' {} {} {}'.format(adapters_file, *reads_filepaths)

        log_filepath = self._log_filepath('fastq-mcf')
        ProgramCaller(command).run(log_filepath=log_filepath)

        return trimmed_reads_filepaths

    def align_to_reference(self, reads_filepaths, ref='GRCh37'):
        """
        Expects a list of two files: forward and reverse reads of the same
        sample. It will search for a reference genome defined by Config.
        """
        params_dict = self.params['bwa']
        params = ['-{} {}'.format(k, v) for k, v in params_dict.items()]
        params_str = ' '.join(params).format(**{
            'reference_genome': Resource('reference_genome')
        })
        command = '{} {} {} {}'.format(self.executables['bwa'], params_str,
                                       *reads_filepaths)
        # redirect stdout to samfile and stderr  to logfile
        log_filepath = self._log_filepath('bwa')
        sam_filepath = join(self.results_dir, self.sample.id + '.sam')
        ProgramCaller(command).run(stdout_sink=sam_filepath,
                                   log_filepath=log_filepath)

    def add_or_replace_read_groups(self, sample):
        params = self.params['AddOrReplaceReadGroups']
        params = ['{}={}'.format(k, v) for k, v in params.items()]
        params_str = ' '.join(params).format(**{
            'sample_id': sample.id,
            'library_id': sample.sequencing.library_id,
            'ngs_id': sample.sequencing.id,
            'input': sample._files('sam'),
            'output': sample._files('bam'),
        })
        command = '{} AddOrReplaceReadGroups {}'.format(
            self.executables['picard-tools'], params_str)
        log_filepath = self._log_filepath('AddOrReplaceReadGroups')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def realign_indels(self, bam_filepath):
        gatk = GATK(bam_filepath)
        gatk.realign_indels()

    def recalibrate_quality_scores(self, bam_filepath):
        gatk = GATK(bam_filepath)
        gatk.recalibrate_quality_scores()

    def call_variants(self, bam_filepath):
        gatk = GATK(bam_filepath)
        gatk.create_vcf()
        gatk.create_gvcf()

    def _log_filepath(self, label):
        return join(self.results_dir, label + '.log')
