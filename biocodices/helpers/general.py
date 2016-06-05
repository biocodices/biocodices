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
