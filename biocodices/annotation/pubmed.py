from functools import lru_cache

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

    @staticmethod
    @lru_cache(maxsize=3000)
    def citation_from_pmid(pmid):
        """
        Fetches a nice AMA citation given a PubMed ID. It uses some site
        I found on the web: https://mickschroeder.com/citation/
        """
        citation_url = 'https://mickschroeder.com/api/citation/cite.php?q={}'
        response = requests.get(citation_url.format(pmid))

        if not response.ok:
            response.raise_for_status()

        return response.json()['html_citation']

