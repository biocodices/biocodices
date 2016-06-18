import luigi

from biocodices.components import Cohort, Sample
from biocodices.programs import trim_adapters, BWA, Picard, GATK


class TrimReads(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self):
        return [luigi.LocalTarget(fn) for fn in self.sample.fastqs]
    def run(self): trim_adapters(self.sample.fastqs)
    def output(self):
        self.sample = Sample(self.sample_dir)
        return [luigi.LocalTarget(fn) for fn in self.sample.trimmed_fastqs]


class AlignReads(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return TrimReads(self.sample_dir)
    def run(self):
        trimmed_fastqs = [target.fn for target in self.requires().output()]
        BWA().align_to_reference(trimmed_fastqs)
    def output(self):
        self.sample = Sample(self.sample_dir)
        return luigi.LocalTarget(self.sample.file('sam'))


class AddOrReplaceReadGroups(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return AlignReads(self.sample_dir)
    def run(self):
        Picard().add_or_replace_read_groups(
            sam_path=self.requires().output().fn,
            sample_id=self.sample.id,
            sample_library_id=self.sample.library_id,
            sequencer_run_id=self.sample.sequencer_run_id,
            out_path=self.output().fn,
        )
    def output(self):
        self.sample = Sample(self.sample_dir)
        return luigi.LocalTarget(self.sample.file('raw.bam'))


class RealignReads(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return AddOrReplaceReadGroups(self.sample_dir)
    def run(self):
        gatk = GATK()
        raw_bam = self.requires().output().fn
        realigned_bam = gatk.realign_reads_around_indels(raw_bam)
        gatk.recalibrate_quality_scores(realigned_bam,
                                        out_path=self.output().fn)
    def output(self):
        self.sample = Sample(self.sample_dir)
        return luigi.LocalTarget(self.sample.file('recalibrated.bam'))


class AlignmentMetrics(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return RealignReads(self.sample_dir)
    def run(self):
        recalibrated_bam = self.requires().output().fn
        Picard().generate_alignment_metrics(recalibrated_bam,
                                            out_path=self.output().fn)
    def output(self):
        self.sample = Sample(self.sample_dir)
        return luigi.LocalTarget(self.sample.file('alignment_metrics.tsv'))


class CoverageMetrics(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return RealignReads(self.sample_dir)
    def run(self):
        recalibrated_bam = self.requires().output().fn
        GATK().create_depth_vcf(recalibrated_bam, out_path=self.output().fn)
    def output(self):
        self.sample = Sample(self.sample_dir)
        return luigi.LocalTarget(self.sample.file('depth_stats.vcf'))


class HaplotypeCall(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return RealignReads(self.sample_dir)
    def run(self):
        recalibrated_bam = self.requires().output().fn
        GATK().create_gvcf(recalibrated_bam, out_path=self.output().fn)
    def output(self):
        self.sample = Sample(self.sample_dir)
        return luigi.LocalTarget(self.sample.file('raw_variants.g.vcf'))


class CohortVariantCalling(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        self.cohort = Cohort(self.base_dir)
        return [HaplotypeCall(sample.dir) for sample in self.cohort.samples]
    def output(self):
        pass

            #  if self.args['--metrics']:
                #  pipeline.add_task(self.cohort.plot_alignment_metrics,
                                #  self.cohort.msg('Plot alignment metrics'))
                #  pipeline.add_task(self.cohort.median_coverages,
                                #  self.cohort.msg('Compute median coverages'))

            #  if self.args['--joint-genotyping']:
                #  gvcf_list = [sample.raw_gvcf for sample in self.cohort.samples]
                #  func = partial(self.gatk.joint_genotyping,
                            #  gvcf_list=gvcf_list,
                            #  output_dir=self.cohort.results_dir)
                #  pipeline.add_task(func, self.cohort.msg('Joint genotyping'))

            #  if self.args['--hard-filtering']:
                #  #  func = partial(self.vcf_munger.hard_filtering,
                            #  #  self.cohort.unfiltered_vcf)
                #  #  pipeline.add_task(func, self.cohort.msg('Hard filtering'))

                #  func = partial(self.gatk.filter_genotypes,
                            #  self.cohort.filtered_vcf)
                #  pipeline.add_task(func, self.cohort.msg('Apply genotype filters'))

                #  func = partial(self.vcf_munger.limit_regions,
                            #  self.cohort.geno_filtered_vcf)
                #  pipeline.add_task(func, self.cohort.msg('Limit VCF regions'))


            #  if self.args['--subset']:
                #  task_group = OrderedDict()

                #  for sample in self.cohort.samples:
                    #  task_label = sample.msg('Subset from multisample VCF')
                    #  task = partial(self.bcftools.subset_samples,
                                #  vcf_path=self.cohort.filtered_vcf,
                                #  sample_ids=[sample.id],
                                #  outfile=sample.filtered_vcf)
                    #  task_group[task_label] = task

                #  pipeline.add_multitask(
                    #  task_group, self.cohort.msg('Subset from multisample VCF'),
                    #  n_processes=n_processes)


if __name__ == '__main__':
    luigi.run()
