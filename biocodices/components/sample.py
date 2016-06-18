import os
import pandas as pd
from termcolor import colored

from biocodices.variant_calling import ReadsMunger, VcfMunger
from biocodices.helpers.language import seconds_to_hms_string
from biocodices.helpers import Config


class Sample:
    reads_format = 'fastq'

    def __init__(self, sample_id, cohort):
        """
        Expects a sample_id and a Cohort object.
        It will look for its fastq files in the cohort.data_dir with
        filenames like: <sample_id>.R1.fastq, <sample_id>.R2.fastq
        (for the forward and reverse reads, respectively).
        """
        self.id = sample_id
        self.cohort = cohort
        self.results_dir = os.path.join(self.cohort.results_dir, self.id)
        self._set_data()

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencer_run_id)

    def analyze_and_trim_reads(self):
        for reads_filepath in self.fastqs:
            self.reads_munger.analyze_reads(reads_filepath)

        self.reads_munger.trim_adapters(self.fastqs)

        for trimmed_filepath in self.trimmed_fastqs:
            self.reads_munger.analyze_reads(trimmed_filepath)

    def alignment_metrics(self):
        self.reads_munger.generate_alignment_metrics(self.recalibrated_bam)
        self.get_median_coverage()

    def get_median_coverage(self):
        self.depth_vcf = self._files('realigned.recalibrated.depth_stats.vcf')
        if not os.path.isfile(self.depth_vcf):
            self.depth_vcf = self.reads_munger.depth_vcf(self.recalibrated_bam)
        depth_stats = VcfMunger.read_depth_stats_vcf(self.depth_vcf)
        return pd.Series(depth_stats).median()

    def read_alignment_metrics(self):
        # This is called by the sample's Cohort to plot everything together.
        fn = self._files('realigned.recalibrated.alignment_metrics.tsv')
        df = pd.read_table(fn, sep='\s+', comment='#')
        df['sample'] = self.id
        return df

    def log(self, extension):
        return os.path.join(self.results_dir,
                            '{}.{}.log'.format(self.id, extension))

    def msg(self, msg):
        prefix = colored('[{}]'.format(self.id), 'cyan')
        return '{} {}'.format(prefix, msg)

    def _log_total_time(self, t1, t2):
        elapsed = seconds_to_hms_string((t2 - t1).seconds)
        with open(self.log('time'), 'w') as logfile:
            msg = 'Variant calling for {} took {}.\n'
            logfile.write(msg.format(self.id, elapsed))

    def _files(self, ext):
        if ext in ['fastq', 'trimmed.fastq']:
            return self._reads_files(ext)

        return '{}.{}'.format(os.path.join(self.results_dir, self.id), ext)

    def _reads_files(self, ext):
        location = self.cohort.data_dir
        if ext == 'trimmed.fastq':
            location = self.results_dir
        read_filepath = os.path.join(location, '{}.{}.{}')
        forward_filepath = read_filepath.format(self.id, 'R1', ext)
        reverse_filepath = read_filepath.format(self.id, 'R2', ext)
        if ext == 'fastq':
            if not os.path.isfile(forward_filepath) or \
               not os.path.isfile(reverse_filepath):
                msg = "I couldn't find BOTH R1 and R2 reads: {}, {}"
                raise OSError(msg.format(forward_filepath, reverse_filepath))

        return forward_filepath, reverse_filepath

    def _set_data(self):
        self.library_id = self._get_library_id()
        self.sequencer_run_id = self._get_seq_run_id()
        self.name = self._get_name()
        self.clinic = self._get_clinic()
        self.long_name = '{} ({}) from {}'.format(
            self.name, self.id, self.clinic)
        self.reads_munger = ReadsMunger(self.results_dir)
        self.vcf_munger = VcfMunger()
        self.fastqs = self._files('fastq')
        self.trimmed_fastqs = self._files('trimmed.fastq')
        self.sam = self._files('sam')
        self.bam = self._files('bam')
        self.realigned_bam = self._files('realigned.bam')
        self.recalibrated_bam = self._files('realigned.recalibrated.bam')
        self.raw_vcf = self._files(Config.filenames['sample_raw_vcf'])
        self.raw_gvcf = self._files(Config.filenames['sample_raw_gvcf'])
        self.filtered_vcf = \
            self._files(Config.filenames['sample_filtered_vcf'])


    # The following methods use two yml files to get stuff that is in the
    # database. Consequently, they have to be updated regularly.
    # It would be better to just read the database, but this is probably
    # faster for the program to start.
    # The sole reason why we need this info (library id, sequencing id,) is
    # for the AddOrRelaceReadGroups step of the variant calling.
    # These methods also mean some hardcoding of biocidces company logic.
    # FIXME: think of a better way of doing this.
    def _get_library_id(self):
        try:
            return Config('db')[self.id][0]
        except KeyError:
            return self.cohort.id

    def _get_seq_run_id(self):
        try:
            return Config('db')[self.id][1]
        except KeyError:
            return self.cohort.id

    def _get_name(self):
        try:
            return Config('db_names')[self.id][0]
        except KeyError:
            return self.id

    def _get_clinic(self):
        try:
            return Config('db_names')[self.id][1]
        except KeyError:
            return 'unknown clinic'
