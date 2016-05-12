# import pandas as pd
from os import makedirs
from os.path import isdir, join, isfile


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
        self.reads_filenames = self._reads_filenames()
        if not isdir(self.results_dir):
            makedirs(self.results_dir, exist_ok=True)

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencing.id)

    def _reads_filenames(self):
        read_filepath = join(self.sequencing.results_dir, '{}.R1.{}')
        forward_filepath = read_filepath.format(self.id, self.__class__.reads_format)
        reverse_filepath = forward_filepath.replace('R1', 'R2')
        if not isfile(forward_filepath) or not isfile(reverse_filepath):
            msg = "I couldn't find both R1 and R2 reads: {}, {}"
            raise OSError(msg.format(forward_filepath, reverse_filepath))

        return forward_filepath, reverse_filepath
