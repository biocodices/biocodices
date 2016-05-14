# import pandas as pd
from os import makedirs
from os.path import isdir, join, isfile

from variant_calling.reads_munger import ReadsMunger


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

    def call_variants(self):
        reads_munger = ReadsMunger(self.results_dir)
        for reads_filepath in self.files['reads']:
            reads_munger.analyze_reads(reads_filepath, self.log('fastqc'))
        trimmed_reads = self.trim_adapters(self.files['reads'],
                                           self.log('fastq-mcf'))
        reads_munger.analyze_reads(trimmed_reads, self.log('trimmed.fastqc'))
        self.files['trimmed_reads'] = trimmed_reads

        #  self.align_to_reference(self.files['trimmed_reads'])

        #  self.variant_call()

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencing.id)

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
