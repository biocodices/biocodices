import luigi
from biocodices.components import Cohort


class TrimReads(luigi.Task):
    cohort_base_dir = luigi.Parameter()

    def run(self):
        self.cohort = Cohort(self.cohort_base_dir)
        for sample in self.cohort.samples:
            sample.analyze_and_trim_reads

    def output(self):
        out = []
        for sample in self.cohort.samples:
            out += [luigi.LocalTarget(reads_file)
                    for reads_file in sample.trimmed_fastqs]
