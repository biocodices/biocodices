from .annotator_with_cache import AnnotatorWithCache
from .my_variant_annotator import MyVariantAnnotator
from .ensembl_parser import EnsemblParser
from .myvariant_parser import MyvariantParser
from .variant_annotator import VariantAnnotator
from .gene_annotation import annotate_gene
from .ensembl import Ensembl
from .dbsnp import DbSNP

from Bio import Entrez
from os.path import isfile

# We need to set an e-mail address for Entrez services
# Warnings about usage will be sent to that address before
# attempting a ban, so it will be nice to receive the warn.
_fn = '~/.mail_address_for_Entrez'
if isfile(_fn):
    with open(_fn) as f:
        Entrez.email = f.read().strip()
else:
    print('WARNING: Please set a mail for Entrez in %s' % _fn)
    print('WARNING: Entrez will notify that mail before banning you')
    print('WARNING: in case your usage is too high.')
    Entrez.email = 'johndoe@example.org'

from .clinvar_variation import ClinvarVariation
from .clinvar_accession import ClinvarAccession
