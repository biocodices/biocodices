from os.path import dirname, join

from biocodices.helpers import Config
from biocodices.programs import AbstractGenomicsProgram, ProgramCaller


class BcfTools(AbstractGenomicsProgram):
    def __init__(self):
        self.label = 'bcftools'
        super(self.__class__, self).__init__(self.label)

    def run(self, module_name, infile, outfile, params_str=None):
        command = '{} {} {} {}'.format(self.executable, module_name,
                                       (params_str or ''), infile)

        log_label = '{}_{}'.format(self.label, module_name)
        log_filepath = join(dirname(outfile), log_label)
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)

    def filter_samples(self, vcf_path, sample_ids, outfile):
        """
        Makes a new VCF keeping the columns you specify.
        - vcf_path: an absolute path to a [g]VCF file.
        - sample_ids: a list of the columns (i.e. samples ids) to keep.
        - outfile: an absolute path for the new vcf file.
        """
        params_str = '--samples ' + ','.join(sample_ids)
        self.run('view', vcf_path, outfile, params_str=params_str)

    def limit_samples(self, gz_vcf_path):
        """
        Subset a VCF keeping the regions contained in the BED defined as
        'panel_amplicons' in the resources.yml file.
        """
        pass

    @staticmethod
    def tabix(gz_vcf_path):
        """Index a bgzip-zipped VCF.gz to let some bcftools modules use it."""
        executable = Config('executables')['tabix']
        command = '{} {}'.format(executable, gz_vcf_path)
        log_filepath = join(dirname(gz_vcf_path), 'tabix')
        ProgramCaller(command).run(log_filepath=log_filepath)
        return gz_vcf_path + '.tbi'

    @staticmethod
    def bgzip(vcf_path):
        """Zips a VCF to vcf.gz to let some bcftools modules use it."""
        executable = Config('executables')['bgzip']
        outfile = vcf_path.replace('.vcf', '.vcf.gz')
        command = '{} -ci {}'.format(executable, vcf_path)
        log_filepath = join(dirname(outfile), 'bgzip')
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
        return outfile
