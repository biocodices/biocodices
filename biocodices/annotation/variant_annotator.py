import requests
import json
import time
import redis
from multiprocessing import Pool
import pandas as pd
from myvariant import MyVariantInfo

from biocodices.annotation import MyvariantParser, EnsemblParser
from biocodices.helpers.general import restful_api_query, in_groups_of


class VariantAnnotator:
    # FIXME: this should check if Redis is present in the system!
    redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

    @classmethod
    def ensembl_cache_get(cls, rs_list):
        ret = {}
        for rs in rs_list:
            info = cls.redis_client.get(rs)
            if info:
                info = json.loads(info.decode('utf8'))
                ret.update({rs: info})
        return ret

    @classmethod
    def ensembl_cache_set(cls, rs_info_dict):
        expire_time = 60 * 60 * 24 * 30  # One month
        for rs, info in rs_info_dict.items():
            if not info: continue
            info = json.dumps(info)
            cls.redis_client.setex(rs, expire_time, info)

    def ensembl_batch_query(self, rs_list, build='GRCh37', use_cache=True):
        """
        Takes a list of rs IDs and queries Ensembl.
        Returns a dictionary with the rs IDs as keys.
        It can annotate with builds GRCh37 and GRCh38.
        """
        rs_list = set([rs for rs in rs_list if rs])
        url = 'http://rest.ensembl.org/variation/homo_sapiens'
        if build == 'GRCh37':
            url = url.replace('rest.', 'grch37.rest.')

        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json'}
        ret = {}
        if use_cache:
            ret.update(self.ensembl_cache_get(rs_list))
            rs_list = rs_list - ret.keys()
            # ^ Remove the annotated ones from cache

        # Ensembl API won't take goups of > 1000 identifiers
        for rs_group in in_groups_of(50, rs_list):
            print('Query Ensembl for %s identifiers...' % len(rs_group))
            payload = json.dumps({'ids': rs_group})
            response = requests.post(url, headers=headers, data=payload)

            if response.ok:
                ret.update(response.json())
                self.ensembl_cache_set(response.json())
                print(' -> 200 OK!')
            else:
                sleep_time = int(response.headers['X-RateLimit-Reset'])
                print(' -> 400 not OK. Sleeping %s seconds...' % sleep_time)
                print('    (Will retry after sleep).')
                time.sleep(sleep_time)
                self.ensembl_batch_query(rs_list)
                # response.raise_for_status()

        return ret

    @classmethod
    def annotate_in_batch(cls, identifiers, processes=10, timeout=10,
                          myvariant=True, ensembl=True):
        """Query MyVariant.info and Ensembl for a list of rs IDs or strings
        like 'chr1:1234:1235'. Runs in parallel. Returns a dictionary with:
            'detail': a merged pandas DataFrame with the SNPs as indices.
            'publications': a list of dicts, each dict is a publication with
                            phenotype, title, and ID.
        """
        with Pool(processes) as pool:
            annotate_args = dict(myvariant=myvariant, ensembl=ensembl)
            results = [pool.apply_async(cls.annotate, (identifier,), annotate_args)
                       for identifier in set(identifiers)]
            results = [result.get(timeout=timeout) for result in results]

        annotations = pd.concat([ann['detail'] for ann in results],
                                ignore_index=True)
        #  index = ['dbsnp_chrom', 'dbsnp_hg19_start', 'rs_id', 'hgvs_id']
        #  annotations.set_index(index, inplace=True)
        publications = {ann['identifier']: ann['publications'] for ann in results}
        return {'detail': annotations, 'publications': publications}

    @classmethod
    def annotate(cls, identifier, myvariant=True, ensembl=True,
                 myvariant_fields=['all']):
        """
        Query MyVariant.info and Ensembl for info about an rs ID. It also
        accepts identifiers like 'chr1:1234-1235'.
        Returns a dict with the following keys:
        'myvariant': a pandas DF with selected info from MyVariant.info.
        'ensembl': a dict with selected info from Ensembl.
        'publications': a list of publications that mention this rs/identifier
                        taken both from GRASP via Myvariant and Ensembl.
        """
        if myvariant:
            myvariant_df, myvariant_publications = \
                cls.parse_myvariant(cls.query_myvariant(identifier))
            if myvariant_df.empty:
                print('No info from MyVariant.info for %s' % identifier)
        if ensembl:
            ensembl_df, ensembl_publications = \
                cls.parse_ensembl(cls.query_ensembl(identifier))

        if myvariant and ensembl:
            annotation = myvariant_df.join(ensembl_df)
            publications = myvariant_publications + ensembl_publications
        elif myvariant:
            annotation = myvariant_df
            publications = myvariant_publications
        elif ensembl:
            annotation = ensembl_df
            publications = ensembl_publications

        annotation['rs_id'] = identifier  # Leaving this for backwards compatibility
        annotation['query'] = identifier
        return {
            'identifier': identifier,
            'detail': annotation,
            'publications': publications
        }

    @staticmethod
    def query_myvariant(identifier, fields=['all']):
        return MyVariantInfo().query(identifier, fields=fields)

    @staticmethod
    def parse_myvariant(results):
        summary = MyvariantParser.parse_query_results(results)
        publications = MyvariantParser.parse_publications(results)
        return summary, publications

    @staticmethod
    def query_ensembl(rs):
        url = "http://rest.ensembl.org"
        url += "/variation/human/{}?phenotypes=1".format(rs)
        try:
            results = restful_api_query(url)
        except requests.exceptions.HTTPError as error:
            # Ensembl API returns an error 400 when the query gives no results.
            print(error)
            results = {}
        return results

    @staticmethod
    def parse_ensembl(results):
        summary = EnsemblParser.parse_query_results(results)
        publications = EnsemblParser.parse_publications(results)
        return summary, publications

    @staticmethod
    def query_cellbase(chromosome, position, allele, key='snp_phenotype'):
        """Queries CellBase restful API. Default key 'snp_phenotype' will fetch
        info about the pheno effect of a SNP. Other options are
        'mutation_phenotype' for COSMIC database and 'consequence_type'."""
        url = 'http://ws.bioinfo.cipf.es/cellbase/rest/latest/hsa/genomic/variant/'
        url += '{chromosome}:{position}:{allele}'
        url += '/{key}?of=json'

        return restful_api_query(url)
