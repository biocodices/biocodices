import requests
import redis
import time
from bs4 import BeautifulSoup
from urllib.error import HTTPError
from multiprocessing import Pool

from biocodices.helpers.general import in_groups_of


class Clinvar:
    # FIXME: this should check if Redis is present in the system!
    # FIXME: redis config should be read from a YML
    _redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

    def annotate(self, cln_ids, use_cache=True, use_web=True, parallel=5,
                 sleep_time=10):
        """
        Annotate a list of Clinvar Variation IDs (also accepts a single ID).
        Not to be confused with Allele IDs!

        When use_cache is True, it will prioritize using Redis cache to get
        info from the given ID. When use_web is True, it will get info from
        Clinvar's API. The priority is on the cache, unless explicitely
        inactivated. It returns a dict where the keys are the passed IDs.

        * parallel: processes to spawn in parallel when querying the web.
        * sleep_time: time to sleep between queries.

        Example:
            > clinvar = Clinvar()
            > clinvar.annotate(['rs12345', 'rs234'])
            # => {'rs12345': { ... }, 'rs234': { ... }}
        """
        if type(cln_ids) == str:
            cln_ids = [cln_ids]

        cln_ids = set([cln_id for cln_id in cln_ids if cln_id])
        xml_dict = {}

        if use_cache:
            for cln_id in cln_ids:
                if self._cache(cln_id):
                    xml_dict[cln_id] = self._cache(cln_id)

            print('Found %s/%s in Clinvar cache' % (len(xml_dict), len(cln_ids)))
            cln_ids = cln_ids - xml_dict.keys()

        if use_web:
            xmls_from_api = self._batch_query(cln_ids, parallel, sleep_time)
            xml_dict.update(xmls_from_api)

        info_dict = {}
        for cln_id, xml in xml_dict.items():
            info_dict.update({cln_id: self._extract_info_from_xml(xml)})

        return info_dict

    @staticmethod
    def _key(cln_id):
        return 'clinvar:%s' % cln_id

    def _url(self, cln_id):
        base_url = 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        params = '?db=clinvar&rettype=variation&id=%s' % cln_id
        return base_url + params

    def _cache(self, cln_id):
        xml_response = self._redis_client.get(self._key(cln_id))
        if not xml_response == b'None':
            return xml_response.decode('utf8')

    def _batch_query(self, cln_ids, parallel, sleep_time):
        """
        Query Clinvar's API for a list of Clinvar Variatino IDs. It returns a
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
                print('Query Clinvar for %s rs IDs' % len(ids_group))
                print(' %s ... %s' % (ids_group[0], ids_group[-1]))
                for cln_id in set(ids_group):
                    results[cln_id] = pool.apply_async(self._query, (cln_id,))

                for cln_id, result in results.items():
                    xml = result.get(timeout=20)
                    if xml:  # Don't save empty dicts
                        xml_dict[cln_id] = xml

        return xml_dict

    def _query(self, cln_id):
        expiration_time = 60 * 60 * 24 * 30  # One month in seconds
        url = self._url(cln_id)
        print('Visiting:', url)

        response = requests.get(url)
        if response.status_code == 200:
            # print(' -> [ 200 OK ] %s...' % response.text[100:120])
            xml_response = response.text
            # print(' Set cache %s' % key(cln_id))
            self._redis_client.setex(self._key(cln_id), expiration_time,
                                     xml_response)
            xml_response = self._redis_client.get(self._key(cln_id))
            return xml_response.decode('utf8')
        else:
            raise HTTPError('Status code was', response.status_code)

    def _extract_info_from_xml(self, raw_xml):
        parsed_xml = raw_xml.split('\n')[1:]
        parsed_xml = '\n'.join(parsed_xml)
        soup = BeautifulSoup(parsed_xml, 'lxml-xml')

        for variation in soup.find_all('VariationReport'):
            variation_id = variation['VariationID']

            rs_list = []
            starts_G37 = []
            stops_G37 = []
            chroms_G37 = []

            for allele in variation.find_all('Allele'):
                rs_elements = allele.find_all('XRef', attrs={'DB': 'dbSNP', 'Type': 'rs'})

                for rs_element in rs_elements:
                    rs_list.append('rs' + rs_element['ID'])

                for seq_location in allele.find_all('SequenceLocation', attrs={'Assembly': 'GRCh37'}):
                    starts_G37.append(seq_location['start'])
                    stops_G37.append(seq_location['stop'])
                    chroms_G37.append(seq_location['Chr'])

            if not starts_G37: print('No starts for', variation_id)
            if not stops_G37: print('No stops for', variation_id)

            rs_id = '|'.join(rs_list) or None # don't assign empty string to rs_id
            start = '|'.join(starts_G37)
            stop = '|'.join(stops_G37)
            chrom = '|'.join(chroms_G37)

            if not rs_id:
                rs_id = '%s:%s-%s' % (chrom, start, stop)
                print('No rs IDs for', variation_id, '-> Assign %s' % rs_id)

            return {'rs_id': rs_id,
                    'chrom_G37': chrom,
                    'start_G37': start,
                    'stop_G37': stop}
