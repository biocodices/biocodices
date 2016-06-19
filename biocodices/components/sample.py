from os.path import join, basename, abspath, expanduser, isfile
from shutil import move
import pandas as pd
from termcolor import colored

from biocodices.variant_calling import VcfMunger
from biocodices.helpers import Config


class Sample:
    reads_format = 'fastq'

    def __init__(self, base_dir):
        """
        Expects a base directory named after a sample. The directory basename
        will be used as a label or ID for the sample. The directory should
        contain fastq files like: <sample_id>.R1.fastq, <sample_id>.R2.fastq
        (for the forward and reverse reads, respectively).
        """
        self.id = basename(base_dir)
        self.dir = abspath(expanduser(base_dir))
        self.fastqs = self._reads_files('fastq')
        self.trimmed_fastqs = self._reads_files('trimmed.fastq')

    def __repr__(self):
        return '<Sample {} from {}>'.format(self.id, self.sequencer_run_id)

    def msg(self, msg):
        prefix = colored('[{}]'.format(self.id), 'cyan')
        return '{} {}'.format(prefix, msg)

    def file(self, ext):
        return '{}.{}'.format(join(self.dir, self.id), ext)

    def _reads_files(self, ext):
        read_filepath = join(self.dir, '{}.{}.{}')
        forward_filepath = read_filepath.format(self.id, 'R1', ext)
        reverse_filepath = read_filepath.format(self.id, 'R2', ext)
        if ext == 'fastq':
            if not isfile(forward_filepath) or \
               not isfile(reverse_filepath):
                msg = "I couldn't find *both* R1 and R2 reads: {}, {}"
                raise FileNotFoundError(msg.format(forward_filepath,
                                                   reverse_filepath))

        return forward_filepath, reverse_filepath

    # The following methods use two yml files to get stuff that is in the
    # database. Consequently, they have to be updated regularly.
    # It would be better to just read the database, but this is probably
    # faster for the program to start.
    # The sole reason why we need this info (library id, sequencing id,) is
    # for the AddOrRelaceReadGroups step of the variant calling.
    # These methods also mean some hardcoding of biocidces company logic.
    # FIXME: think of a better way of doing this.
    @property
    def library_id(self):
        try:
            return Config('db')[self.id][0]
        except KeyError:
            return self.id

    @property
    def sequencer_run_id(self):
        try:
            return Config('db')[self.id][1]
        except KeyError:
            return self.id

    @property
    def name(self):
        try:
            return Config('db_names')[self.id][0]
        except KeyError:
            return self.id

    @property
    def clinic(self):
        try:
            return Config('db_names')[self.id][1]
        except KeyError:
            return 'unknown clinic'
