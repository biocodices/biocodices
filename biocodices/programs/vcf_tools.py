from os.path import dirname, join
from biocodices.programs import AbstractGenomicsProgram, ProgramCaller
from biocodices.helpers import Config


class VcfTools(AbstractGenomicsProgram):
    def __init__(self):
        super(self.__class__, self).__init__('vcftools')

    @staticmethod
    def subset(vcf_path, sample_ids, outfile):
        """
        Makes a new VCF keeping the columns you specify.
        - vcf_path: an absolute path to a [g]VCF file.
        - sample_ids: a list of the columns (i.e. samples ids) to keep.
        - outfile: an absolute path for the new vcf file.
        """
        command = Config.executables['vcf-subset']
        for sample_id in sample_ids:
            command += ' -c {}'.format(sample_id)
        command += ' {}'.format(vcf_path)

        log_filepath = join(dirname(outfile), 'vcf-subset.log')
        ProgramCaller(command).run(stdout_sink=outfile,
                                   log_filepath=log_filepath)
