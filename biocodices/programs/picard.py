from os.path import dirname, join
from biocodices.helpers import Config, Resource
from biocodices.programs import ProgramCaller
from biocodices.helpers.general import rename_tempfile


class Picard:
    def __init__(self):
        self.executable = Config('executables')['picard-tools']
        self.params = Config('parameters')['picard-tools']
        self.reference_genome = Resource('reference_genome')
        self.known_variants = Resource('dbsnp:GRCh37')

    def run(self, module_name, infile, outfile, extra_params_variables={},
            extra_output_extension=None):
        params_dict = self.params[module_name]
        params = ['{}={}'.format(k, v) for k, v in params_dict.items()]
        params_variables = {
            'reference_genome': self.reference_genome,
            'known_variants': self.known_variants,
            'input': infile,
            'output': outfile + '.temp',
        }
        params_variables.update(extra_params_variables)
        params_str = ' '.join(params).format(**params_variables)
        command = '{} {} {}'.format(self.executable, module_name, params_str)
        log_filepath = join(dirname(infile), module_name)
        ProgramCaller(command).run(log_filepath=log_filepath)
        rename_tempfile(outfile, extra_output_extension)

    def alignment_metrics(self, infile):
        outfile = infile.replace('.bam', '.alignment_metrics.tsv')
        self.run('CollectAlignmentSummaryMetrics', infile, outfile)

    def variant_calling_metrics(self, infile):
        outfile = infile.replace('.vcf', '.variant_calling_metrics.tsv')
        self.run('CollectVariantCallingMetrics', infile, outfile)

    def add_or_replace_read_groups(self, sample):
        extra_params_variables = {
            'sample_id': sample.id,
            'library_id': sample.library_id,
            'ngs_id': sample.sequencer_run_id,
        }
        self.run('AddOrReplaceReadGroups', sample.sam, sample.bam,
                 extra_params_variables, 'bai')
