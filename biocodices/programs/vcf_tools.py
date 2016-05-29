from os import rename
from os.path import dirname, join
from biocodices.programs import AbstractGenomicsProgram, ProgramCaller
from biocodices.helpers.general import rename_tempfile


class VcfTools(AbstractGenomicsProgram):
    def __init__(self):
        super(self.__class__, self).__init__('vcftools')
