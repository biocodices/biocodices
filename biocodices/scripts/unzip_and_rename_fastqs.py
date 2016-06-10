#!/usr/bin/env python

import re
import subprocess
import sys
from os.path import join, expanduser, isdir, isfile, abspath
from os import rename, makedirs
from glob import glob
from shutil import copy2


def main(data_dir):
    # Pattern to extract the sample ID
    filename_pattern = r'(SAR\d+).*01-(\d{5}|Control).*(R\d)'

    gzipped_fastq_glob = join(data_dir, '*.fastq.gz')
    gzipped_fastq_filepaths = glob(gzipped_fastq_glob)

    msg = 'Backuping the {} zipped fastqs and renaming the copies...'
    print(msg.format(len(gzipped_fastq_filepaths)))

    for fp in gzipped_fastq_filepaths:
        # Copy original files to archive directory
        original_files_dir = join(data_dir, 'zipped_fastqs')
        makedirs(original_files_dir, exist_ok=True)
        copy2(fp, original_files_dir)

    # This is a separate loop so we make sure the backuping is complete
    # before attempting the unzipping.
    for fp in gzipped_fastq_filepaths:
        # Rename with a cleaner sample ID
        match = re.search(filename_pattern, fp)
        inta_sample_id, sample_id, read = match.groups(1)
        if 'Control' in fp:
            new_fp = join(data_dir, 'Control_{}.{}.fastq.gz'.format(inta_sample_id, read))
        else:
            new_fp = join(data_dir, 'ENPv1_{}.{}.fastq.gz'.format(sample_id, read))
            if isfile(new_fp):
                msg = ('File {} already exists. Repeated sample ID? Check the '
                       'input fastq files for sample {}/{}.')
                raise Exception(msg.format(new_fp, inta_sample_id, sample_id))
        rename(fp, new_fp)

    print('Unzipping!')
    for fp in glob(join(data_dir, "*.fastq.gz")):  # Get the renamed filepaths
        command = 'gzip -d {}'.format(fp)
        print(command)
        subprocess.run(command.split(' '), check=True)

    print('\n=> Done. Go check {}'.format(abspath(data_dir)))


if __name__ == '__main__':
    args = sys.argv
    if len(args) != 2:
        print('\nUsage:')
        print('{} <data dir with fastq.gz files>'.format(__file__))
        sys.exit()

    data_dir = expanduser(sys.argv[1])
    if not isdir(data_dir):
        print("\nAre you sure that directory exists?\n{}".format(data_dir))

    main(data_dir)
