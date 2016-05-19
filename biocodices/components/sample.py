import logging
from datetime import datetime
from os import makedirs
from os.path import isdir, join, isfile

from biocodices.variant_calling.reads_munger import ReadsMunger


class Sample:
    reads_format = 'fastq'

    def __init__(self, sample_id, sequencing):
        """
        Expects a sample_id and a Sequencing object.
        It will look for its fastq files in the sequencing.data_dir with
        filenames like: <sample_id>.R1.fastq, <sample_id>.R2.fastq
        (for the forward and reverse reads, respectively).
        """
        self.id = sample_id
        self.sequencing = sequencing
        self.results_dir = join(self.sequencing.results_dir, self.id)

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencing.id)

    def call_variants(self):
        if not isdir(self.results_dir):
            makedirs(self.results_dir, exist_ok=True)

        t1 = datetime.now()
        munger = ReadsMunger(self.id, self.results_dir)

        #  for reads_filepath in self._files('fastq'):
            #  munger.analyze_reads(reads_filepath)

        munger.trim_adapters(self._files('fastq'))
        #  for trimmed_filepath in self._files('trimmed.fastq'):
            #  munger.analyze_reads(trimmed_filepath)

        # TODO: implement this
        # munger.multiqc(self._files('fastq') + self._files('trimmed.fastq'))

        munger.align_to_reference(self._files('trimmed.fastq'))
        munger.add_or_replace_read_groups(self)
        # TODO: delete the bamfile

        #  self.variant_call()
        t2 = datetime.now()
        self._log_total_time(t1, t2)

    def log(self, extension):
        return join(self.results_dir, '{}.{}.log'.format(self.id, extension))

    def _log_total_time(self, t1, t2):
        timedelta = (t2 - t1).seconds
        with open(self.log('time'), 'w') as logfile:
            msg = 'Variant calling for {} took {} seconds.\n'
            logfile.write(msg.format(self.id, timedelta))

    def _files(self, ext):
        if ext in ['fastq', 'trimmed.fastq']:
            setattr(self, ext, self._reads_files(ext))
            return getattr(self, ext)

        filepath = '{}.{}'.format(join(self.results_dir, self.id), ext)
        setattr(self, ext, filepath)
        return getattr(self, ext)

    def _reads_files(self, ext):
        location = self.sequencing.data_dir
        if ext == 'trimmed.fastq':
            location = self.results_dir
        read_filepath = join(location, '{}.{}.{}')
        forward_filepath = read_filepath.format(self.id, 'R1', ext)
        reverse_filepath = read_filepath.format(self.id, 'R2', ext)
        if not isfile(forward_filepath) or not isfile(reverse_filepath):
            msg = "I couldn't find BOTH R1 and R2 reads: {}, {}"
            raise OSError(msg.format(forward_filepath, reverse_filepath))

        return forward_filepath, reverse_filepath
