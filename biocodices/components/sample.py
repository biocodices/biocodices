from os import walk, getcwd
from os.path import join, basename, isfile


class Sample:
    reads_format = 'fastq'

    def __init__(self, sample_id):
        """
        Expects a sample_id that will *exactly* match the name of ONLY ONE
        subdirectory under the current working dir.
        The ideal use case is that you run bioco from the root dir of a Cohort,
        so then this class will just search for its subdirectory in the results
        folder.
        """
        self.id = sample_id
        self.dir = self._find_subdir(self.id)
        self.fastqs = self._reads_files('fastq')
        self.trimmed_fastqs = self._reads_files('trimmed.fastq')

    def file(self, ext):
        return '{}.{}'.format(join(self.dir, self.id), ext)

    def __repr__(self):
        return '<Sample {}>'.format(self.id)

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

    @staticmethod
    def _find_subdir(sample_id):
        matching_dirs = []
        for dirpath, _, _ in walk(getcwd()):
            if sample_id == basename(dirpath):
                matching_dirs.append(dirpath)

        if len(matching_dirs) > 1:
            msg = 'More than one subdirectory has this ID: %s' % sample_id
            raise ValueError(msg)
        elif len(matching_dirs) == 0:
            msg = 'No subdirectory found with this ID: %s' % sample_id
            raise ValueError(msg)


        return matching_dirs[0]
