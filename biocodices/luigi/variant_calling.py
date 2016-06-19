"""
     _     _                     _ _
    | |__ (_) ___   ___ ___   __| (_) ___ ___  ___
    | '_ \| |/ _ \ / __/ _ \ / _` | |/ __/ _ \/ __|
    | |_) | | (_) | (_| (_) | (_| | | (_|  __/\__ \\
    |_.__/|_|\___/ \___\___/ \__,_|_|\___\___||___/


Usage:
    bioco TASK_NAME [options]
    bioco (-h | --help)

Options:
    -d --database DB_NAME               Database name to use for annotation
                                        of the variatns in the task
                                        DatabaseAnnotation. If you run this
                                        task without providing this option,
                                        the process will exit with an error.
"""

from docopt import docopt
import luigi

from biocodices import software_name
from biocodices.helpers import logo
from biocodices.components import Cohort, Sample
from biocodices.programs import trim_adapters, BWA, Picard, GATK
from biocodices.variant_calling import VcfMunger


class UnzipAndCoppyFastqs(luigi.ExternalTask):
    sample_dir = luigi.Parameter()
    def output(self):
        self.sample = Sample(self.sample_dir)
        return [luigi.LocalTarget(fn) for fn in self.sample.fastqs]


class TrimReads(luigi.Task):
    sample_dir = luigi.Parameter()
    def requires(self): return UnzipAndCoppyFastqs(self.sample_dir)
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
        vcf_munger = VcfMunger()
        raw_bam = self.requires().output().fn
        realigned_bam = vcf_munger.realign_reads_around_indels(raw_bam)
        vcf_munger.recalibrate_quality_scores(realigned_bam,
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


class CohortAlignmentMetrics(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return [AlignmentMetrics(sample.dir) for sample in self.cohort.samples]
    def run(self):
        self.cohort.plot_alignment_metrics(
            metrics_files=[task.output().fn for task in self.requires()],
            out_path=self.output().fn)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        return luigi.LocalTarget(self.cohort.file('alignment_metrics.png'))


class JointGenotyping(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return [HaplotypeCall(sample.dir) for sample in self.cohort.samples]
    def run(self):
        gvcf_list = [task.output().fn for task in self.requires()]
        GATK().joint_genotyping(gvcf_list=gvcf_list, out_path=self.output().fn)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        return luigi.LocalTarget(self.cohort.file('raw_variants.vcf'))


class HardFiltering(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return JointGenotyping(self.cohort.dir)
    def run(self):
        raw_vcf = self.requires().output().fn
        VcfMunger().hard_filtering(raw_vcf, out_path=self.output().fn)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        return luigi.LocalTarget(self.cohort.file('filtered.vcf'))


class GenotypeFiltering(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return HardFiltering(self.cohort.dir)
    def run(self):
        raw_vcf = self.requires().output().fn
        GATK().filter_genotypes(raw_vcf, out_path=self.outfile)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        self.infile = self.requires().output().fn
        self.outfile = self.infile.replace('.vcf', '.geno.vcf')
        return luigi.LocalTarget(self.outfile)


class LimitRegions(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return GenotypeFiltering(self.cohort.dir)
    def run(self):
        VcfMunger().limit_regions(self.infile, out_path=self.outfile)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        self.infile = self.requires().output().fn
        self.outfile = self.infile.replace('.vcf', '.lim.vcf')
        return luigi.LocalTarget(self.outfile)


class SnpEffAnnotation(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return LimitRegions(self.cohort.dir)
    def run(self):
        VcfMunger().annotate_with_snpeff(self.infile, out_path=self.outfile)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        self.infile = self.requires().output().fn
        self.outfile = self.infile.replace('.vcf', '.Eff.vcf')
        return luigi.LocalTarget(self.outfile)


class VEPAnnotation(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    def requires(self):
        return SnpEffAnnotation(self.cohort.dir)
    def run(self):
        VcfMunger().annotate_with_VEP(self.infile, out_path=self.outfile)
    def output(self):
        self.cohort = Cohort(self.base_dir)
        self.infile = self.requires().output().fn
        self.outfile = self.infile.replace('.vcf', '.VEP.vcf')
        return luigi.LocalTarget(self.outfile)


class DatabaseAnnotation(luigi.Task):
    base_dir = luigi.Parameter(default='.')
    database = luigi.Parameter()
    def requires(self):
        return VEPAnnotation(self.cohort.dir)
    def run(self):
        VcfMunger().annotate_pubmed_with_DB(
            vcf_path=self.infile,
            database=self.database,
            citations_table='variationscitations',
            rs_column='VariationName',
            pubmed_column='PubmedID',
            out_path=self.outfile,
        )
    def output(self):
        self.cohort = Cohort(self.base_dir)
        self.infile = self.requires().output().fn
        self.outfile = self.infile.replace('.vcf', '.db.vcf')
        return luigi.LocalTarget(self.outfile)


def run_pipeline():
    docopt(__doc__, version=software_name)
    print(logo())
    welcome_msg = 'Welcome to {}! Starting the Luigi pipeline...\n'
    print(welcome_msg.format(software_name))
    luigi.run()

if __name__ == '__main__':
    run_pipeline()
