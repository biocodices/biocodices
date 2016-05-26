from shutil import move
from os import makedirs
from os.path import join, exists, dirname, abspath
from itertools import product


def params_dict_to_str(params_dict):
    params = ['-{} {}'.format(k, v) for k, v in params_dict.items()]
    return ' '.join(params)


def rename_tempfile(outfile, extra_file_extension=None):
    move(outfile + '.temp', outfile)
    if extra_file_extension:
        extra_file = outfile + '.temp.{}'.format(extra_file_extension)
        move(extra_file, extra_file.replace('.temp', ''))


def touch_all_the_logs(base_dir, dir_list):
    # I wrote this just to be able to run a `tail -f *.log` on every log
    # during the variant calling, even for logs that don't yet exist but that
    # will be created during the process. It's a necessarily hardcoded list,
    # I guess:
    [makedirs(d, exist_ok=True) for d in dir_list]
    path = abspath(join(dirname(__file__), '../scripts/log_filenames.txt'))
    with open(path, 'r') as f:
        log_filenames = [l.strip() for l in f.readlines()]
    log_filepaths = [join(d, fn + '.log')
                     for d, fn in product(dir_list, log_filenames)]
    # log_filepaths += [join(base_dir, 'GenotypeGVCFs.log')]
    for log_filepath in log_filepaths:
        if not exists(log_filepath):
            open(log_filepath, 'a').close()


def biocodices_logo():
    path = abspath(join(dirname(__file__), '../scripts/logo.txt'))
    with open(path, 'r') as logo_file:
        logo_string = logo_file.read()
    return logo_string
