import time
from multiprocessing import Pool

from Bio import Entrez
from biocodices.annotation import AnnotatorWithCache
from biocodices.helpers import in_groups_of


class ClinvarAccession(AnnotatorWithCache):
    @staticmethod
    def _key(rs):
        return 'clinvar:accession:%s' % rs

    @staticmethod
    def _query(accession_number):
        search_term = '%s[ClinVar accession]' % accession_number
        handle = Entrez.esearch(db='clinvar', term=search_term)
        response = Entrez.read(handle)
        return response['IdList']

    def _batch_query(self, accn_list, parallel, sleep_time):
        info_dict = {}

        with Pool(parallel) as pool:
            for i, accn_group in enumerate(in_groups_of(parallel, accn_list)):
                group_results = {}

                for accession_number in accn_group:
                    group_results[accession_number] = pool.apply_async(
                            self._query, (accession_number, ))

                for accession_number, result in group_results.items():
                    id_list = result.get()
                    info_dict[accession_number] = id_list

        self._cache_set(info_dict)
        return info_dict

