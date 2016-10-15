from Bio import Entrez
from bs4 import BeautifulSoup
from multiprocessing import Pool

from biocodices.annotation import AnnotatorWithCache


class ClinvarVariation(AnnotatorWithCache):
    @staticmethod
    def _key(cln_id):
        return 'clinvar:variation:%s' % cln_id

    @staticmethod
    def _url(cln_id):
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        params = '?db=clinvar&rettype=variation&id=%s' % cln_id
        return base_url + params

    @staticmethod
    def _query(cln_id):
        handle = Entrez.efetch(db='clinvar', id=cln_id, rettype='variation')
        return handle.read()

    @classmethod
    def parse_xml_dict(cls, xml_dict):
        """
        Parses a dict like { 'RCV1234': '<xml ...>' }. Meant for the result
        of #annotate() for this class.
        Returns a dict like { 'RCV1234': { data about the accession } }
        """
        info_dict = {}
        with Pool(6) as pool:
            results = {key: pool.apply_async(cls.parse_xml, (xml, ))
                       for key, xml in xml_dict.items()}

            for key, result in results.items():
                info_dict[key] = result.get()

        return info_dict

    @staticmethod
    def parse_xml(xml):
        """
        Parses a dict like { '671234': '<xml ...>' }. Meant for the result
        of #annotate() for this class. Returns a dict like:
        { '671234': { data about the accession } }
        """
        parsed_xml = '\n'.join(xml.split('\n')[1:])
        soup = BeautifulSoup(parsed_xml, 'lxml-xml')

        variations = soup.find_all('VariationReport')
        if len(variations) > 1:
            raise Exception('Many variations in this XML')

        variation = variations[0]
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
            for pubmed in clinical_significance.find_all('ID', attrs={'Source': 'PubMed'}):
                pubmed_ids.add(pubmed.text)

        ret = {
            'chrom_G37': '|'.join(chroms_G37),
            'start_G37': '|'.join(starts_G37),
            'stop_G37': '|'.join(stops_G37),
            'pubmed_ids': '|'.join(list(pubmed_ids)),
            'variation_id': variation_id
        }

        rs_id = '|'.join(rs_list)
        if not rs_id:
            rs_id = '%s:%s-%s' % (ret['chrom_G37'], ret['start_G37'], ret['stop_G37'])
        ret['rs_id'] = rs_id

        return ret
