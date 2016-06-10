#!/usr/bin/env python

"""
 _     _                     _ _
| |__ (_) ___   ___ ___   __| (_) ___ ___  ___
| '_ \| |/ _ \ / __/ _ \ / _` | |/ __/ _ \/ __|
| |_) | | (_) | (_| (_) | (_| | | (_|  __/\__ \\
|_.__/|_|\___/ \___\___/ \__,_|_|\___\___||___/

Usage:
    bioco -d BASE_DIR [options]
    bioco (-h | --help)

Options:
    -d --base-dir BASE_DIR              Base directory for the cohort.
                                        Should contain a "data" subdirectory
                                        with .fastq R1 and R2 files named
                                        after each sample.
    -c --complete-pipeline              Do the whole pipeline.
                                        Equivalent to options: -tavjfm.
    -n --processes <N>                  Number of max simultaneous processes
                                        to spawn in some parts of the pipeline
                                        that have been parallelized.
    -t --trim-reads                     Trim adapters from the fastqs.
    -a --align-reads                    Align reads to the reference genome.
    -v --create-vcf                     Call variants for each sample (from
                                        the .bam to a .g.vcf).
    -j --joint-genotyping               Do a joint genotyping for the cohort.
                                        From samples gVCFs to a single VCF.
    -f --hard-filtering                 Do a hard filtering on the cohort VCF.
                                        Parameters will be read from:
                                        ~/.biocodices/parameters.yml
    -s --samples SAMPLE SAMPLE...       Do the steps *only* on these samples.
    -k --skip-samples SAMPLE SAMPLE...  Skip these samples.
    -m --metrics                        Generate and plot alignment and
                                        coverage metrics.
    -h --help                           Show this screen.

Examples:
    bioco (-d DIR) [--complete-pipeline] [-n N_PROCESSES]

"""

from docopt import docopt

from biocodices import software_name
from biocodices.components import PipelineCreator


def cli():
    """Runs the pipeline according to the passed CLI arguments."""
    arguments = docopt(__doc__, version=software_name)

    if arguments['--processes']:
        arguments['--processes'] = int(arguments['--processes'])
    else:
        arguments['--processes'] = 1

    if arguments['--complete-pipeline']:
        arguments['--trim-reads'] = True
        arguments['--align-reads'] = True
        arguments['--create-vcf'] = True
        arguments['--metrics'] = True
        arguments['--joint-genotyping'] = True
        arguments['--hard-filtering'] = True

    pipeline_creator = PipelineCreator(arguments)
    pipeline_creator.pre_pipeline()
    pipeline = pipeline_creator.build_pipeline()
    pipeline.run()
    pipeline_creator.post_pipeline()

if __name__ == '__main__':
    cli()
