import subprocess
from os.path import dirname, join, isfile


class Plink:
    def __init__(self, bfile_path):
        if not isfile(bfile_path + '.bed'):
            raise FileNotFoundError(bfile_path + '.bed')
        self.input_bfile = bfile_path  # .replace('.bed', '')
        self.workdir = dirname(bfile_path)

    def __repr__(self):
        return '<Plink for "{}">'.format(self.input_bfile)

    def ped(self):
        return self.run('--recode', make_bed=False) + '.ped'

    def extract(self, snps_filename, out):
        return self.run('--extract {}'.format(snps_filename), out=out)

    def keep_fam(self, famfile, out):
        return self.run('--keep-fam {}'.format(famfile), out=out)

    def fst(self, clusters_file, out=None):
        return self.run('--within {} --fst'.format(clusters_file), out=out,
                        make_bed=False)

    def assoc(self, adjust=True):
        options = '--assoc'
        if adjust:
            options += ' --adjust'
        return self.run(options, make_bed=False)

    def run(self, options, out=None, make_bed=True):
        if out is None:
            out = self.input_bfile
        out = join(self.workdir, out)
        command_template = 'plink --bfile {} {}'
        if make_bed:
            command_template += ' --make-bed'
        command_template += ' --out {}'
        command = command_template.format(self.input_bfile, options, out)
        subprocess.run(command.split(' '))

        return out

    @staticmethod
    def bim_fields():
        return ["chr", "rs_id", "morgans", "position", "A1", "A2"]

    @staticmethod
    def fam_fields():
        return ['family', 'sample', 'father', 'mother', 'sexcode', 'phenotype']

    @staticmethod
    def make_bed_from_ped(path_label):
        command_template = 'plink --file {} --make-bed --out {}'
        command = command_template.format(path_label, path_label)
        subprocess.run(command.split(' '))

        return path_label

