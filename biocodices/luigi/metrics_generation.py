import luigi

from biocodices.components import Cohort, Sample
from biocodices.programs import Picard, GATK
from biocodices.luigi import RealignReads


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


class CohortAlignmentMetrics(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        for sample in self.cohort.samples: yield AlignmentMetrics(sample.dir)
    def run(self):
        self.cohort.plot_alignment_metrics(
            metrics_files=[task.output().fn for task in self.requires()],
            out_path=self.output().fn)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        return luigi.LocalTarget(self.cohort.file('alignment_metrics.png'))


class CohortCoverageMetrics(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        for sample in self.cohort.samples: yield CoverageMetrics(sample.dir)
    def run(self): pass  # TODO: implement
    def output(self): pass  # TODO: implement
