from os.path import dirname, join
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
