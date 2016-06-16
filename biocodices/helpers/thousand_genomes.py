from tempfile import NamedTemporaryFile

from biocodices.helpers import Resource
from biocodices.components import Dataset
from biocodices.programs import Plink


def dataset_from_1000Genomes_at_loci(rs_list, out_label):
    """Extracts a list of markers from the 1000 Genomes VCFs.
    Makes a new set of plink bed/bim/fam files with the resulting
    dataset (that set of markers for all the 1000G samples),
    named after the given 'label-path'."""

    # Create a file with the rs_list, one per line
    with NamedTemporaryFile() as tempfile:
        for rs in rs_list: tempfile.write(rs + '\n')

        # Run plink extract with that file
        #  Plink.extract_snps_from_vcf(snps_file=tempfile.name,
                                    #  vcf_path=)

    # Read the out files and create a dataset
    # Use the 1000 Genomes populations as index?


