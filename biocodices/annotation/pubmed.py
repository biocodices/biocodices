import json
from Bio import Entrez

from biocodices.annotation import AnnotatorWithCache


class PubMed(AnnotatorWithCache):
    @staticmethod
    def _key(pmid):
        return 'pubmed:{0}'.format(pmid)

    @staticmethod
    def _query(pmid):
        handle = Entrez.efetch(id=pmid, db='pubmed', rettype='xml')
        response = Entrez.read(handle)
        # The response is a Bio.Entrez.Parser.ListElement that cannot be
        # pickled, so it raises errors with the multiprocessing.
        # We need to jsonify it.
        return json.dumps(response)

