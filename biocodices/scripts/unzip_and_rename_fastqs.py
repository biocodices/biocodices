#!/usr/bin/env python

import re
import subprocess
import sys
from os import rename, makedirs
from os.path import basename, dirname, join, isfile, abspath
from glob import glob
from shutil import move


# Pattern to extract the sample ID from the fastq filenames
PATTERNS = [
    r'(SAR\d+).*01-(\d{5}).*_(R\d)',
    r'(SAR\d+).*01(\d{5}).*_(R\d)',
    r'(sar\d+).*ENPv1_(\d{5}).*_(R\d)',
]

# Add this prefix to the new FASTQ filenames
PREFIX_TO_ADD = 'ENPv1'


def main():
    print('\nUnzipping the fastq.gz files\n')
    for fn in glob('*.fastq.gz'):
        if any(bit in fn for bit in ['cont', 'Cont']):
            continue  # Skip controls

        command = 'gzip -dc {}'.format(fn)
        new_fn = fn.replace('.fastq.gz', '.fastq')

        with open(new_fn, 'w') as new_file:
            print(command + ' > ' + new_fn)
            subprocess.run(command.split(' '), stdout=new_file, check=True)

    print('\nRenaming the unzipped fastqs with sample IDs\n')
    for fn in glob('*.fastq'):
        matched_this_one = False

        for pattern in PATTERNS:
            match = re.search(pattern, fn)
            if not match:
                continue

            inta_sample_id, sample_id, read_number = match.groups(1)
            new_fn = '{}_{}.{}.fastq'.format(PREFIX_TO_ADD, sample_id,
                                             read_number)
            if isfile(new_fn):
                msg = ('(!) File {} already exists.\n\n'
                       'Repeated sample ID? Check '
                       'the input fastq files for sample {}/{}.')
                print(msg.format(new_fn, inta_sample_id, sample_id))
                sys.exit()

            print(fn + ' -> ' + new_fn)
            rename(fn, new_fn)

            matched_this_one = True
            break  # if a pattern was found, don't keep trying with more

        if not matched_this_one:
            print("(!) Please add a PATTERN to match this {}\n".format(fn))
            sys.exit()

    #  print('\nMoving the renamed files to a sister directory "results"\n')
    #  for fn in glob('*.fastq'):
        #  current_filepath = abspath(fn)
        #  dest_dir = join(dirname(current_filepath), '../results')
        #  dest_dir = abspath(dest_dir)
        #  makedirs(dest_dir, exist_ok=True)
        #  move(current_filepath, join(dest_dir, basename(current_filepath)))

    #  print("=> Done. You can find the unzipped fastq files in ../restuls\n")


if __name__ == '__main__':
    main()
