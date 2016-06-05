from functools import partial


#  class Component(object):
    #  def __init__(self, *args, **kwargs):
        #  raise NotImplementedError('Abstract class, do not instantiate.')

    #  def run(self):
        #  raise NotImplementedError


#  class Task(Component):
    #  def __init__(self, *args, **kwargs):
        #  Component.__init__(self, *args, **kwargs)

    #  def run(self):
        #  raise NotImplementedError


# class Pipeline(Component):
class Pipeline:
    def __init__(self, *args, **kwargs):
        # Component.__init__(self, *args, **kwargs)
        self.tasks = []

    def add_task(self, task):
        self.tasks.append(task)

    def remove_task(self, task):
        self.tasks.remove(task)

    def run(self):
        #  last_output = None
        for task in self.tasks:
            output = task(last_output)
            #  last_output = output


class PipelinePrepare:
    pipeline = Pipeline()

    def run(self, cohort,
            trim_reads=True, align_reads=True,
            create_vcfs=True, joint_genotyping=True,
            hard_filtering=True, plot_metrics=False):

        if trim_reads or align_reads or create_vcfs:
            for sample in cohort.samples:
                task = partial(sample.call_variants, trim_reads=trim_reads,
                               align_reads=align_reads, create_vcfs=create_vcfs)
                # split this assignment in mini-sample-tasks too!!!!
                pipeline.add_task(task)

        if plot_metrics:
            # self.printlog('Plot some alignment metrics for the cohort.')
            pipeline.add_task(cohort.plot_alignment_metrics)
            # self.printlog('Compute median coverage of the cohort.')
            pipeline.add_task(cohort.median_coverages)

        if joint_genotyping:
            # self.printlog('Joint genotyping.')
            pipeline.add_task(cohort.joint_genotyping)

        if hard_filtering:
            # self.printlog('Hard filtering the multisample VCF')
            cohort.apply_filters_to_vcf(cohort.unfiltered_vcf)

            # self.printlog('Split the multisample VCF per sample')
            for sample in cohort.samples:
                self.vcf_munger.filter_samples(cohort.filtered_vcf, [sample.id],
                                               sample.filtered_vcf)
