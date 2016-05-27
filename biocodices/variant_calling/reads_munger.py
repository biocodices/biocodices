from os.path import join, basename

from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller, GATK, Picard, BWA
from biocodices.helpers.general import params_dict_to_str


class ReadsMunger:
    def __init__(self, sample, results_dir):
        self.sample = sample
        self.results_dir = results_dir
        self.executables = Config('executables')
        self.params = Config('parameters')
        self.picard = Picard()
        self.gatk = GATK()
        self.bwa = BWA()

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__, self.sample.id)

    def analyze_reads(self, reads_filepath):
        command = '{} {} -o {}'.format(self.executables['fastqc'],
                                       reads_filepath, self.results_dir)
        log_filepath = self._file('fastqc')
        ProgramCaller(command).run(log_filepath=log_filepath)

    def trim_adapters(self, reads_filepaths):
        """
        Expects a list of two files: forward and reverse reads of the same
        sample. It will search for an adapters file defined by Config.
        """
        trimmed_fastqs = []
        for reads_filepath in reads_filepaths:
            new_fn = basename(reads_filepath).replace('.fastq', '.trimmed.fastq')
            new_filepath = join(self.results_dir, new_fn)
            trimmed_fastqs.append(new_filepath)

        params_str = params_dict_to_str(self.params['fastq-mcf'])
        for trimmed_fastq in trimmed_fastqs:
            params_str += ' -o {}'.format(trimmed_fastq)

        command = '{} {}'.format(self.executables['fastq-mcf'], params_str)
        adapters_file = Resource('illumina_adapters_file')
        command += ' {} {} {}'.format(adapters_file, *reads_filepaths)

        log_filepath = self._file('fastq-mcf')
        ProgramCaller(command).run(log_filepath=log_filepath)

        return trimmed_fastqs

    def align_to_reference(self, reads_filepaths):
        # NOTE: the human referece genome used is defined in the config yml
        # ~/.biocodices/resources.yml as "&reference_genome_default".
        outfile = self.bwa.align_to_reference(reads_filepaths)
        return outfile

    def add_or_replace_read_groups(self, sample):
        self.picard.add_or_replace_read_groups(sample)

    def realign_reads_around_indels(self, bam_filepath):
        self.gatk.realign_reads_around_indels(bam_filepath)

    def recalibrate_quality_scores(self, realigned_bam):
        self.gatk.recalibrate_quality_scores(realigned_bam)

    def alignment_metrics(self, recalibrated_bam):
        self.picard.alignment_metrics(recalibrated_bam)

    def create_vcf(self, recalibrated_bam):
        self.gatk.create_vcf(recalibrated_bam)

    def create_gvcf(self, recalibrated_bam):
        self.gatk.create_gvcf(recalibrated_bam)

    def _file(self, label):
        return join(self.results_dir, label)
