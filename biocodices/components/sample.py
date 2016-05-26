from os import makedirs, remove
from os.path import isdir, join, isfile
from datetime import datetime
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
        self.results_dir = join(self.cohort.results_dir, self.id)
        self._set_data()

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencer_run_id)

    def call_variants(self, trim_reads=True, align_reads=True,
                      create_vcfs=True):
        if not isdir(self.results_dir):
            makedirs(self.results_dir, exist_ok=True)

        t1 = datetime.now()
        print(colored('* {}'.format(self.long_name), 'yellow'))
        print('Result files and logs will be put in:')
        print(self.results_dir, '\n')

        #  if trim_reads:
            #  self.analyze_and_trim_reads()
        #  if align_reads:
            #  self.align_reads()
            #  self.process_alignment_files()
        if create_vcfs:
            self.create_variant_files()

        t2 = datetime.now()
        self._log_total_time(t1, t2)
        print()

    def analyze_and_trim_reads(self):
        self.printlog('Analyze reads')
        for reads_filepath in self.fastqs:
            self.reads_munger.analyze_reads(reads_filepath)

        self.printlog('Trim adapters')
        self.reads_munger.trim_adapters(self.fastqs)

        self.printlog('Analyze trimmed reads')
        for trimmed_filepath in self.trimmed_fastqs:
            self.reads_munger.analyze_reads(trimmed_filepath)

    def align_reads(self):
        self.printlog('Align reads to reference')
        self.reads_munger.align_to_reference(self.trimmed_fastqs)

    def process_alignment_files(self):
        #  self.printlog('Add or replace read groups')
        #  self.reads_munger.add_or_replace_read_groups(self)

        #  self.printlog('Delete sam file')
        #  remove(self._files('sam'))

        #  self.printlog('Realign indels')
        #  self.reads_munger.realign_indels(self.bam)

        #  self.printlog('Recalibrate read quality scores')
        #  self.reads_munger.recalibrate_quality_scores(self.bam)

        self.printlog('Generate alignment metrics')
        self.reads_munger.alignment_metrics(self.recalibrated_bam)

    def create_variant_files(self, vcf=True, gvcf=True):
        #  if vcf:
            #  self.printlog('Create vcf')
            #  self.reads_munger.create_vcf(self.bam)
        if gvcf:
            self.printlog('Create gvcf')
            self.reads_munger.create_gvcf(self.bam)

    def apply_filters_to_vcf(self):
        if isfile(self.joint_vcf):
            input_vcf = self.joint_vcf
        else:
            msg = ("WARNING: I couldn't find a vcf from joint genotyping, "
                   "so I'll use the regular one: {}".format(self.vcf))
            self.printlog(msg)
            input_vcf = self.vcf

        variant_files = self.vcf_munger.create_snp_and_indel_vcfs(input_vcf)
        self.snps_vcf, self.indels_vcf = variant_files

        self.printlog('Apply SNP filters')
        self.snps_vcf = self.vcf_munger.apply_filters(self.snps_vcf, 'snps')

        self.printlog('Apply indel filters')
        self.indels_vcf = self.vcf_munger.apply_filters(self.indels_vcf,
                                                        'indels')
        self.printlog('Merge the filtered vcfs')
        self.vcf_munger.merge_variant_vcfs([self.snps_vcf, self.indels_vcf],
                                            outfile=self.filtered_vcf)

    def log(self, extension):
        return join(self.results_dir, '{}.{}.log'.format(self.id, extension))

    def printlog(self, msg):
        timestamp = datetime.now().strftime('%H:%M:%S')
        prefix = colored('[{}][{}]'.format(timestamp, self.id), 'blue')
        print('{} {}'.format(prefix, msg))

    def _log_total_time(self, t1, t2):
        elapsed = seconds_to_hms_string((t2 - t1).seconds)
        with open(self.log('time'), 'w') as logfile:
            msg = 'Variant calling for {} took {}.\n'
            logfile.write(msg.format(self.id, elapsed))

    def _files(self, ext):
        if ext in ['fastq', 'trimmed.fastq']:
            return self._reads_files(ext)

        return '{}.{}'.format(join(self.results_dir, self.id), ext)

    def _reads_files(self, ext):
        location = self.cohort.data_dir
        if ext == 'trimmed.fastq':
            location = self.results_dir
        read_filepath = join(location, '{}.{}.{}')
        forward_filepath = read_filepath.format(self.id, 'R1', ext)
        reverse_filepath = read_filepath.format(self.id, 'R2', ext)
        if ext == 'fastq':
            if not isfile(forward_filepath) or not isfile(reverse_filepath):
                msg = "I couldn't find BOTH R1 and R2 reads: {}, {}"
                raise OSError(msg.format(forward_filepath, reverse_filepath))

        return forward_filepath, reverse_filepath

    def _get_library_id(self):
        return Config('db')[self.id][0]

    def _get_seq_run_id(self):
        return Config('db')[self.id][1]

    def _get_name(self):
        return Config('db_names')[self.id][0]

    def _get_clinic(self):
        return Config('db_names')[self.id][1]

    def _set_data(self):
        self.library_id = self._get_library_id()
        self.sequencer_run_id = self._get_seq_run_id()
        self.name = self._get_name()
        self.clinic = self._get_clinic()
        self.long_name = '{} ({}) from {}'.format(
            self.name, self.id, self.clinic)
        self.reads_munger = ReadsMunger(self, self.results_dir)
        self.vcf_munger = VcfMunger(self, self.results_dir)
        self.fastqs = self._files('fastq')
        self.trimmed_fastqs = self._files('trimmed.fastq')
        self.sam = self._files('sam')
        self.bam = self._files('bam')
        self.recalibrated_bam = self._files('realigned.recalibrated.bam')
        self.vcf = self._files('vcf')
        self.gvcf = self._files('g.vcf')
        self.joint_vcf = self._files('joint.vcf')
        self.filtered_vcf = self._files('filtered.vcf')
