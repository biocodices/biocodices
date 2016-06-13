import re
import yaml
from os.path import dirname, join, isfile, basename

from biocodices.helpers import Config
from biocodices.programs import ProgramCaller


class Plink:
    bim_fields = ["chr", "rs_id", "morgans", "position", "A1", "A2"]
    fam_fields = ['FID', 'IID', 'father', 'mother', 'sexcode', 'phenotype']
    executable = Config('executables')['plink']

    def __init__(self, bfile_path):
        self.tests = self.config('association_tests')

        self.label_path = bfile_path
        self.bed = self.label_path + '.bed'
        self.bim = self.label_path + '.bim'
        self.fam = self.label_path + '.fam'

        for ext in ['.bed', '.bim', '.fam']:
            fn = self.label_path + ext
            if not isfile(fn):
                msg = "I couldn't find the {} file: {}"
                raise FileNotFoundError(msg.format(ext, fn))

        self.workdir = dirname(self.label_path)

    def __repr__(self):
        return '<{} for "{}">'.format(self.__class__.__name__,
                                      basename(self.label_path))

    def make_ped(self):
        return self.run('--recode')

    def make_traw(self):
        return self.run('--recode A-transpose')

    def extract(self, snps_filename, out):
        """Keeps the variants with IDs listed in the 'snps_filename' passed.
        Creates a new set of plink bed/bim/fam files the chosen 'out' label."""
        return self.run('--extract {}'.format(snps_filename), out=out)

    def keep_fam(self, famfile, out):
        """Keeps samples with the family IDs listed in the passed famfile.
        Generates a new set of plink bed/bim/fam files."""
        return self.run('--keep-fam {}'.format(famfile), out=out)

    def fst(self, clusters_file, out=None):
        return self.run('--within {} --fst'.format(clusters_file), out=out)

    def pre_tests_filter(self, mind=0.1, maf=0.05, geno=0.1, hwe=0.0001,
                         out=None):
        """Run a conventional set of filters previous to tests. Produces a new
        set of plink bed/bim/fam files on which the tests should be ran."""
        params = '--mind {} --maf {} --geno {} --hwe {}'
        params = params.format(mind, maf, geno, hwe)
        return self.run(params, out=out, make_bed=True)

    def association_test(self, test_name, famfile=None, out=None):
        test = self.tests[test_name]
        if famfile:
            self.fam = famfile
            out = out or famfile.replace('.fam', '')
        return self.run(test['params'], out=out)

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
        """Run raw CLI options on the bed/bim/fam set in self. Optionally, make
        a new bed file with 'out' label."""
        command = '{} --bed {} --bim {} --fam {}'
        # ^ WARNING: Keep the trailing space
        command = command.format(self.executable, self.bed, self.bim, self.fam)
        command += ' {}'.format(options)
        if make_bed:
            command += ' --make-bed'
        out_path = join(self.workdir, (out or self.label_path))
        command += ' --out {}'.format(out_path)
        ProgramCaller(command).run()
        return self._out_filepath_from_log(out_path)

    def _out_filepath_from_log(self, out_label):
        with open(out_label + '.log', 'r') as logfile:
            # Take the last lines from the logfile:
            lines = [line.strip() for line in logfile.readlines()][-4:]
        pattern = r'({}.*?)(\s|$)'.format(out_label)
        matches = [re.search(pattern, line) for line in lines]
        # Assuming a single match for the filepath in the last lines of the log
        return [match for match in matches if match][0].group(0).strip()

    @classmethod
    def make_bed_from_ped(cls, path_label):
        """Make a plink bed file from a plink ped file."""
        command_template = 'plink --file {} --make-bed --out {}'
        command = command_template.format(path_label, path_label)
        ProgramCaller(command).run()

        return path_label

    @classmethod
    def make_bed_from_filtered_vcf(cls, path_label, out_label=None):
        """Makes a bed file from a vcf, leaving out the markers that didn't
        pass previous filteres (i.e. don't have PASS in the filter column)."""
        command = '{} --vcf {} --vcf-filter --make-bed'
        command = command.format(cls.executable, path_label)
        out_label = out_label or basename(path_label).replace('.vcf', '')
        command += ' --out {}'.format(out_label)
        ProgramCaller(command).run()
        return out_label

    @classmethod
    def extract_snps_from_vcf(cls, snps_file, vcf_path, out_label):
        """Extract the SNPs listed in a file from a source VCF. Keeps all the
        samples. Writes the bed/bim/fam files to out_label."""
        command = '{} --vcf {} --extract {} --make-bed --out {}'
        command = command.format(cls.executable, vcf_path, snps_file,
                                 out_label)
        # out_label = out_label or basename(vcf_path).replace('.vcf', '')
        command += ' --out {}'.format(out_label)
        ProgramCaller(command).run()
        return out_label

    @staticmethod
    def config(label=None):
        config_file = join(dirname(__file__), 'plink.yml')
        with open(config_file, 'r') as fh:
            config_dict = yaml.load(fh)
        if label is None:
            return config_dict

        return config_dict[label]
