from datetime import datetime
from shutil import move
from itertools import product
from os.path import join
import requests
import sys


def touch_all_the_logs(cohort):
    """
    Creates empty files for all the logs that will be written during the
    variant calling, so you can tail -f them.
    """
    # I wrote this just to be able to run a `tail -f *.log` on every log
    # during the variant calling, even for logs that don't yet exist but
    # that would later be created. What I do is just creating the empty
    # logs beforehand. It's a necessarily hardcoded list, I guess:
    samples_dirs = [sample.dir for sample in cohort.samples]

    for path, fn in product(samples_dirs, sample_log_filenames()):
        log_path = join(path, fn + '.log')
        open(log_path, 'w').close()

    for fn in cohort_log_filenames():
        log_path = join(cohort.results_dir, fn + '.log')
        open(log_path, 'w').close()


def timestamp(sep=':', hour=True, date=False):
    template = '%H{0}%M{0}%S'.format(sep)
    if date:
        template = '%Y-%m-%d_{}'.format(template)
    return datetime.now().strftime(template)


def params_dict_to_str(params_dict):
    params = ['-{} {}'.format(k, v) for k, v in params_dict.items()]
    return ' '.join(params)


def rename_tempfile(outfile, extra_file_extension=None):
    move(outfile + '.temp', outfile)
    if extra_file_extension:
        extra_file = outfile + '.temp.{}'.format(extra_file_extension)
        move(extra_file, extra_file.replace('.temp', ''))


def restful_api_query(url):
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    return response.json()


def sample_log_filenames():
    return [
        'fastqc',
        'fastq-mcf',
        'bwa',
        'AddOrReplaceReadGroups',
        'RealignerTargetCreator',
        'IndelRealigner',
        'BaseRecalibrator',
        'PrintReads',
        'CollectAlignmentSummaryMetrics',
        'DiagnoseTargets',
        'HaplotypeCaller_gvcf',
        'bcftools_view_samples',
    ]


def cohort_log_filenames():
    return [
        'GenotypeGVCFs',
        'SelectVariants_indels',
        'SelectVariants_snps',
        'VariantFiltration_indels_filter',
        'VariantFiltration_snps_filter',
        'CombineVariants',
        'VariantFiltration_genotype_filter',
        'SnpEff',
        'VEP',
    ]


def logo():
    return """
   _     _                     _ _
  | |__ (_) ___   ___ ___   __| (_) ___ ___  ___
  | '_ \| |/ _ \ / __/ _ \ / _` | |/ __/ _ \/ __|
  | |_) | | (_) | (_| (_) | (_| | | (_|  __/\__ \\
  |_.__/|_|\___/ \___\___/ \__,_|_|\___\___||___/

"""
