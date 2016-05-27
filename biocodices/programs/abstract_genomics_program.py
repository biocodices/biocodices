from biocodices.helpers import Config, Resource


class AbstractGenomicsProgram:
    def __init__(self, program_name):
        self.executable = Config('executables')[program_name]
        self.params = Config('parameters')[program_name]
        self.reference_genome = Resource('reference_genome')
        self.known_indels = [Resource('indels:1000G'),
                             Resource('indels:mills')]
        self.panel_amplicons = Resource('panel_amplicons')
        self.known_variants = Resource('dbsnp:GRCh37')
