from os.path import join, basename

from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller, GATK
from biocodices.helpers.general import params_dict_to_str, rename_tempfile


class ReadsMunger:
    def __init__(self, sample, results_dir):
        self.sample = sample
        self.results_dir = results_dir
        self.executables = Config('executables')
        self.params = Config('parameters')
        self.gatk = GATK()

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

    def align_to_reference(self, reads_filepaths, ref='GRCh37'):
        """
        Expects a list of two files: forward and reverse reads of the same
        sample. It will search for a reference genome defined by Config.
        """
        outfile = join(self.results_dir, self.sample.id + '.sam')

        params_str = params_dict_to_str(self.params['bwa']).format(**{
            'reference_genome': Resource('reference_genome')
        })
        for reads_file in reads_filepaths:
            params_str += ' {}'.format(reads_file)

        command = '{} {}'.format(self.executables['bwa'], params_str)
        # redirect stdout to samfile and stderr  to logfile
        log_filepath = self._file('bwa')
        ProgramCaller(command).run(stdout_sink=outfile + '.temp',
                                   log_filepath=log_filepath)
        rename_tempfile(outfile)

    def add_or_replace_read_groups(self, sample):
        outfile = sample._files('bam')
        params_dict = self.params['AddOrReplaceReadGroups']
        params = ['{}={}'.format(k, v) for k, v in params_dict.items()]
        params_str = ' '.join(params).format(**{
            'sample_id': sample.id,
            'library_id': sample.library_id,
            'ngs_id': sample.sequencer_run_id,
            'input': sample._files('sam'),
            'output': outfile + '.temp',
        })
        command = '{} AddOrReplaceReadGroups {}'.format(
            self.executables['picard-tools'], params_str)
        log_filepath = self._file('AddOrReplaceReadGroups')
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, 'bai')

    def realign_indels(self, bam_filepath):
        self.gatk.set_bamfile(bam_filepath)  # TODO: improve this
        self.gatk.realign_indels()

    def recalibrate_quality_scores(self, bam_filepath):
        self.gatk.set_bamfile(bam_filepath)  # TODO: improve this
        self.gatk.recalibrate_quality_scores()

    def create_vcf(self, bam_filepath):
        self.gatk.set_bamfile(bam_filepath)  # TODO: improve this
        self.gatk.create_vcf()

    def create_gvcf(self, bam_filepath):
        self.gatk.set_bamfile(bam_filepath)  # TODO: improve this
        self.gatk.create_gvcf()

    def _file(self, label):
        return join(self.results_dir, label)
