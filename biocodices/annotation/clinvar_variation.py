import re
import time
from Bio import Entrez
from bs4 import BeautifulSoup
from multiprocessing import Pool

from biocodices.annotation import AnnotatorWithCache
from biocodices.helpers.general import in_groups_of


class ClinvarVariation(AnnotatorWithCache):
    @staticmethod
    def _key(cln_id):
        return 'clinvar:variation:%s' % cln_id

    @staticmethod
    def _url(cln_id):
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        params = '?db=clinvar&rettype=variation&id=%s' % cln_id
        return base_url + params

    def _batch_query(self, cln_ids, parallel, sleep_time):
        """
        Query Clinvar's API for a list of Clinvar Variation IDs. It returns a
        dict with mappings and some data. Runs in <parallel> processes and
        sleeps <sleep> seconds between queries.
        """
        with Pool(parallel) as pool:
            xml_dict = {}

            for i, ids_group in enumerate(in_groups_of(parallel, cln_ids)):
                if i > 0:
                    print(' Sleep %s seconds' % sleep_time)
                    time.sleep(sleep_time)

                results = {}
                print('Query Clinvar for %s Variation IDs' % len(ids_group))
                print(' %s ... %s' % (ids_group[0], ids_group[-1]))
                for cln_id in set(ids_group):
                    results[cln_id] = pool.apply_async(self._query, (cln_id,))

                for cln_id, result in results.items():
                    xml = result.get(timeout=20)
                    if xml:  # Don't save empty dicts
                        xml_dict[cln_id] = xml

        self._cache_set(xml_dict)
        return xml_dict

    @staticmethod
    def _query(cln_id):
        handle = Entrez.efetch(db='clinvar', id=cln_id, rettype='variation')
        return handle.read()

    @staticmethod
    def parse_xml(xml):
        parsed_xml = '\n'.join(xml.split('\n')[1:])
        soup = BeautifulSoup(parsed_xml, 'lxml-xml')

        variations = soup.find_all('VariationReport')
        if len(variations) > 1:
            raise Exception('Many variations in this XML')

        variation =  variations[0]
        variation_id = variation['VariationID']

        rs_list = []
        starts_G37 = []
        stops_G37 = []
        chroms_G37 = []
        pubmed_ids = set()

        for allele in variation.find_all('Allele'):
            rs_elements = allele.find_all('XRef', attrs={'DB': 'dbSNP', 'Type': 'rs'})

            for rs_element in rs_elements:
                rs_list.append('rs' + rs_element['ID'])

            for seq_location in allele.find_all('SequenceLocation', attrs={'Assembly': 'GRCh37'}):
                starts_G37.append(seq_location['start'])
                stops_G37.append(seq_location['stop'])
                chroms_G37.append(seq_location['Chr'])

        # Publications are listed for the Variation, not for each Allele
        for clinical_significance in variation.find_all('ClinicalSignificance'):
            desc = clinical_significance.find('Description').text

            for pubmed in clinical_significance.find_all('ID', attrs={'Source': 'PubMed'}):
                pubmed_ids.add(pubmed.text)

        ret = {
            'chrom_G37': '|'.join(chroms_G37),
            'start_G37': '|'.join(starts_G37),
            'stop_G37': '|'.join(stops_G37),
            'pubmed_ids': '|'.join(list(pubmed_ids))
        }

        rs_id = '|'.join(rs_list)
        if not rs_id:
            rs_id = '%s:%s-%s' % (ret['chrom'], ret['start'], ret['stop'])
        ret['rs_id'] = rs_id

        return ret
