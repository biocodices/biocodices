import re
from glob import glob
from os.path import expanduser, basename, join, abspath

from .sample import Sample
from biocodices.helpers.language import plural


class SequencerRun:
    def __init__(self, directory):
        """
        A SequencerRun object expects the path to a directory that will have
        'data' and 'results' subdirectories. fastq files from samples will be
        looked for in the 'data' subdir.
        """
        self.dir = abspath(expanduser(directory))
        self.id = basename(self.dir)
        self.data_dir = join(self.dir, 'data')
        self.results_dir = join(self.dir, 'results')
        self.library_id = 'ENPv1_LIB_00001'  ## FIXME: ! Check this!!!!

    def __repr__(self):
        tmpl = '<SequencerRun({})>'
        return tmpl.format(self.dir)

    def __str__(self):
        n_samples = len(self.samples())
        tmpl = 'Sequencer run "{}" with {}'
        return tmpl.format(self.id, plural('sample', n_samples))

    def samples(self):
        glob_expr = join(self.data_dir, '*.{}'.format(Sample.reads_format))
        all_reads_filenames = sorted(glob(glob_expr))
        sample_ids = [re.search(r'(.*).R1', basename(reads_fn)).groups(1)[0]
                      for reads_fn in all_reads_filenames
                      if 'R1' in reads_fn]

        return [Sample(sample_id, self) for sample_id in sample_ids]
