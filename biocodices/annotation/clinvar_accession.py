import re
import time
from multiprocessing import Pool

from bs4 import BeautifulSoup
from Bio import Entrez
from biocodices.annotation import AnnotatorWithCache
from biocodices.helpers import in_groups_of


class ClinvarAccession(AnnotatorWithCache):
    @staticmethod
    def _key(rs):
        return 'clinvar:accession:%s' % rs

    @staticmethod
    def _query(accession_number):
        handle = Entrez.efetch(db='clinvar', rettype='clinvarset',
                               id=accession_number)
        return handle.read()

    def _batch_query(self, accn_list, parallel, sleep_time):
        info_dict = {}

        with Pool(parallel) as pool:
            for i, accn_group in enumerate(in_groups_of(parallel, accn_list)):
                if i > 0:
                    time.sleep(sleep_time)

                print('Query Entrez for %s clinvar accessions' % len(accn_group))
                print(' %s ... %s' % (accn_group[0], accn_group[-1]))
                group_results = {}

                for accession_number in accn_group:
                    group_results[accession_number] = pool.apply_async(
                            self._query, (accession_number, ))

                for accession_number, result in group_results.items():
                    xml = result.get()
                    info_dict[accession_number] = xml
                    # Update the cache as you go along, not at the end
                    self._cache_set({accession_number: xml})

        return info_dict

    @classmethod
    def parse_xml_dict(cls, xml_dict):
        """
        Parses a dict like { 'RCV1234': '<xml ...>' }. Meant for the result
        of #annotate() for this class.
        Returns a dict like { 'RCV1234': { data about the accession } }
        """
        return {key: cls.parse_xml(xml) for key, xml in xml_dict.items()}

    @staticmethod
    def parse_xml(xml):
        """
        Parses an accession XML response from Entrez/Clinvar.
        """
        soup = BeautifulSoup(xml, 'xml')

        variation_ids = []
        var_filter = {'Type': 'Variant', 'ID': re.compile('\d+')}
        for variant in soup.find_all(attrs=var_filter):
            variation_ids.append(variant['ID'])

        omim_ids = []
        mim_filter = {'Type': 'MIM', 'DB': 'OMIM', 'ID': re.compile('\d+')}
        for mim in soup.find_all(attrs=mim_filter):
            omim_ids.append(mim['ID'])

        rs_ids = []
        dbsnp_filter = {'Type': 'rs', 'DB': 'dbSNP', 'ID': re.compile('\d+')}
        for rs in soup.find_all(attrs=dbsnp_filter):
            rs_ids.append('rs' + rs['ID'])

        info = {
            'variation_ids': variation_ids,
            'dbsnp_ids': rs_ids,
            'omim_ids': omim_ids
        }
        return info
