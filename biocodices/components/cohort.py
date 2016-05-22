from biocodices.programs import GATK
from biocodices.helpers.language import plural
from biocodices.variant_calling import VcfMunger


class Cohort:
    def __init__(self, samples):
        """
        Expects a list of Sample objects, not necessarily from the same
        sequencer run.
        """
        self.samples = samples
        self.sequencer_runs = [sample.sequencer_run for sample in self.samples]
        if len(self.samples) == 0:
            raise Exception('Are you trying to create empty Cohort?')

    def __repr__(self):
        return '<{}>'.format(self.__str__())

    def __str__(self):
        n_sequencer_runs = len(self.sequencer_runs)
        n_samples = len(self.samples)
        tmpl = '{} with {} from {}'
        return tmpl.format(self.__class__.__name__,
                           plural('sample', n_samples),
                           plural('sequencer run', n_sequencer_runs))

    def call_variants(self):
        for sample in self.samples:
            sample.call_variants()

        self.joint_vcf = self.joint_genotyping()

        for sample in self.samples:
            sample.joint_vcf = self.joint_vcf
            sample.apply_filters_to_vcf()

    def joint_genotyping(self):
        gatk = GATK()
        gvcf_list = [sample.gvcf for sample in self.samples]
        # Since the samples might come from multiple sequencer runs, we
        # pick on of them and put the joint vcf in its results directory.
        # If they're samples from only one sequence run, the path will be
        # the expected one (the results directory parent to all sample
        # directories).
        output_dir = self.sequencer_runs[0].results_dir
        self.joint_vcf = gatk.joint_genotyping(gvcf_list, output_dir)
        self.split_joint_vcf()

    def split_joint_vcf(self):
        for sample in self.samples:
            outfile = sample._files('post-joint.vcf')
            VcfMunger.subset(vcf=self.joint_vcf, sample_ids=[sample.id],
                             outfile=outfile)
