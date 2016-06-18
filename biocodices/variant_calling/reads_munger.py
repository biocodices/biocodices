from os.path import join, dirname

from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller, GATK, Picard, BWA, FastQC
from biocodices.helpers.general import params_dict_to_str


class ReadsMunger:
    def __init__(self):
        self.fastqc = FastQC()
        self.bwa = BWA()
        self.picard = Picard()
        self.gatk = GATK()

    def trim_adapters(self, reads_filepaths):
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

    def align_to_reference(self, reads_filepaths):
        # NOTE: the human referece genome used is defined in the config yml
        # ~/.biocodices/resources.yml as "&reference_genome_default".
        outfile = self.bwa.align_to_reference(reads_filepaths)
        return outfile

    def process_alignment_files(self, sam_path, sample_id,
                                sample_library_id, sequencer_run_id):
        picard_args = sam_path, sample_id, sample_library_id, sequencer_run_id
        bam = self.picard.add_or_replace_read_groups(*picard_args)
        realigned_bam = self.gatk.realign_reads_around_indels(bam)
        recalibrated_bam = self.gatk.recalibrate_quality_scores(realigned_bam)

        return recalibrated_bam

    def generate_alignment_metrics(self, recalibrated_bam):
        return self.picard.alignment_metrics(recalibrated_bam)

    def depth_vcf(self, recalibrated_bam):
        return self.gatk.create_depth_vcf(recalibrated_bam)

    #  def generate_variant_calling_metrics(self, filtered_vcf):
        #  self.picard.variant_calling_metrics(filtered_vcf)

    def _file(self, label):
        return join(self.results_dir, label)
