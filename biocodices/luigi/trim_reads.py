import luigi

from biocodices.components import Cohort
from biocodices.programs import FastQC
from biocodices.variant_calling import ReadsMunger


class TrimReads(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    reads_munger = ReadsMunger()
    fastqc = FastQC()

    def output(self):
        self.cohort = Cohort(self.base_dir)

        return [luigi.LocalTarget(trimmed_fastq)
                for sample in self.cohort.samples
                for trimmed_fastq in sample.trimmed_fastqs]

    def run(self):
        for sample in self.cohort.samples:

            # The analysis is done for each reads file separatedly
            #  for reads_filepath in sample.fastqs:
                #  self.fastqc.analyze_reads(reads_filepath)

            # The trimming is done for the pair in the same command:
            self.reads_munger.trim_adapters(sample.fastqs)

if __name__ == '__main__':
    luigi.run()
