# import pandas as pd
from datetime import datetime
from os import makedirs
from os.path import isdir, join, isfile

from biocodices.variant_calling.reads_munger import ReadsMunger


class Sample:
    reads_format = 'fastq'

    def __init__(self, sample_id, sequencing):
        """
        This class expects a sample_id and a Sequencing object.
        It will look for its fastq files in the sequencing.data_dir with
        filenames like: <sample_id>.R1.fastq, <sample_id>.R2.fastq
        (for the forward and reverse reads, respectively).
        """
        self.id = sample_id
        self.sequencing = sequencing
        self.results_dir = join(self.sequencing.results_dir, self.id)
        self.files = {}
        self.files['reads'] = self._reads_filepaths()
        if not isdir(self.results_dir):
            makedirs(self.results_dir, exist_ok=True)

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencing.id)

    def call_variants(self):
        t1 = datetime.now()
        munger = ReadsMunger(self.id, self.results_dir)
        for reads_filepath in self.files['reads']:
            munger.analyze_reads(reads_filepath)
        self.files['trimmed_reads'] = munger.trim_adapters(self.files['reads'])
        for trimmed_filepath in self.files['trimmed_reads']:
            munger.analyze_reads(trimmed_filepath)
        munger.align_to_reference(self.files['trimmed_reads'])

        #  self.variant_call()
        t2 = datetime.now()
        self._log_total_time(t1, t2)

    def _reads_filepaths(self):
        read_filepath = join(self.sequencing.data_dir, '{}.R1.{}')
        forward_filepath = read_filepath.format(self.id,
                                                self.__class__.reads_format)
        reverse_filepath = forward_filepath.replace('R1', 'R2')
        if not isfile(forward_filepath) or not isfile(reverse_filepath):
            msg = "I couldn't find both R1 and R2 reads: {}, {}"
            raise OSError(msg.format(forward_filepath, reverse_filepath))

        return forward_filepath, reverse_filepath

    def log(self, extension):
        return join(self.results_dir, '{}.{}.log'.format(self.id, extension))

    def _log_total_time(self, t1, t2):
        timedelta = (t2 - t1).seconds
        with open(self.log('time'), 'w') as logfile:
            msg = 'Variant calling took {} seconds.\n'.format(timedelta)
            logfile.write(msg)

