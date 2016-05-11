import subprocess
import re
from os.path import dirname, join, isfile


class Plink:
    def __init__(self, bfile_path):
        if not isfile(bfile_path + '.bed'):
            raise FileNotFoundError(bfile_path + '.bed')
        self.input_bfile = bfile_path  # .replace('.bed', '')
        self.workdir = dirname(bfile_path)

    def __repr__(self):
        return '<Plink for "{}">'.format(self.input_bfile)

    def make_ped(self):
        return self.run('--recode')

    def make_traw(self):
        return self.run('--recode A-transpose')

    def extract(self, snps_filename, out, make_bed_True):
        return self.run('--extract {}'.format(snps_filename), out=out)

    def keep_fam(self, famfile, out, make_bed=True):
        return self.run('--keep-fam {}'.format(famfile), out=out)

    def fst(self, clusters_file, out=None):
        return self.run('--within {} --fst'.format(clusters_file), out=out)

    def assoc(self, adjust=True):
        options = '--assoc'
        if adjust:
            options += ' --adjust'
        return self.run(options)

    def model(self, mperm=None, fisher=False):
        options = '--model'
        if fisher:
            options += ' --fisher'
        if mperm:
            options += ' --mperm {}'.format(mperm)
        return self.run(options)

    def test_missing(self):
        return self.run('--test-missing')

    def run(self, options, out=None, make_bed=False):
        command_template = 'plink --bfile {} {}'
        command_template += ' --out {}'
        out = out or self.input_bfile
        out = join(self.workdir, out)
        if make_bed:
            command_template += ' --make-bed'
        command = command_template.format(self.input_bfile, options, out)
        self.__class__.execute(command)
        return self._out_filepath_from_log(out)

    def _out_filepath_from_log(self, out_label):
        with open(out_label + '.log', 'r') as logfile:
            # Take the last lines from the logfile:
            lines = [line.strip() for line in logfile.readlines()][-4:]
        pattern = r'({}.*)\b'.format(out_label)
        matches = [re.search(pattern, line) for line in lines]
        # Assuming a single match for the filepath in the last lines of the log
        return [match for match in matches if match][0].group(0)

    @staticmethod
    def execute(command):
        try:
            subprocess.run(command.split(' '), check=True)
        except subprocess.CalledProcessError as error:
            print('Problems with this command:\n\n', ' '.join(error.cmd))
            raise error

    @staticmethod
    def bim_fields():
        return ["chr", "rs_id", "morgans", "position", "A1", "A2"]

    @staticmethod
    def fam_fields():
        return ['family', 'sample', 'father', 'mother', 'sexcode', 'phenotype']

    @classmethod
    def make_bed_from_ped(cls, path_label):
        command_template = 'plink --file {} --make-bed --out {}'
        command = command_template.format(path_label, path_label)
        cls.execute(command)

        return path_label
