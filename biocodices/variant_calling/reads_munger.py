from os.path import join, basename

from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller, GATK, Picard, BWA, FastQC
from biocodices.helpers.general import params_dict_to_str


class ReadsMunger:
    def __init__(self, results_dir):
        self.results_dir = results_dir
        # ^ This class needs a results_dir arg since it can't infer the results
        # directory from the location of the input fastq files, which are in
        # the **data** dir of each project.
        self.executables = Config('executables')
        self.picard = Picard()
        self.gatk = GATK()
        self.bwa = BWA()
        self.fastqc = FastQC()

    def analyze_reads(self, reads_filepath):
        self.fastqc.analyze_reads(reads_filepath, self.results_dir)

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

        params_str = params_dict_to_str(Config.params['fastq-mcf'])
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
