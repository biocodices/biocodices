import time
import json
from multiprocessing import Pool
from Bio import Entrez
from io import StringIO

from biocodices.annotation import AnnotatorWithCache
from biocodices.helpers import in_groups_of


class PubMed(AnnotatorWithCache):
    @staticmethod
    def _key(pmid):
        return 'pubmed:{0}'.format(pmid)

    @staticmethod
    def _query(pmid):
        handle = Entrez.efetch(id=pmid, db='pubmed', rettype='xml')
        # The response is a Bio.Entrez.Parser.ListElement that cannot be
        # pickled, so it raises errors with the multiprocessing.
        # Instead, return the XML
        return handle.read()

    def _batch_query(self, pmids, parallel, sleep_time):
        info_dict = {}

        with Pool(parallel) as pool:
            for pmid_group in in_groups_of(parallel, pmids):

                print('Query Entrez for %s PubMed IDs' % len(pmid_group))
                print(' %s ... %s' % (pmid_group[0], pmid_group[-1]))
                group_results = {pmid: pool.apply_async(self._query, (pmid, ))
                                 for pmid in pmid_group}

                for pmid, result in group_results.items():
                    info = Entrez.read(StringIO(result.get()))
                    info_dict[pmid] = info
                    # Update the cache as you go along, not at the end
                    self._cache_set({pmid: json.dumps(info)})

                if len(pmids) > parallel:
                    time.sleep(sleep_time)

        return info_dict

