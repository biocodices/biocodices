import time
import requests
import json
import redis
from multiprocessing import Pool, TimeoutError
from itertools import chain

from biocodices.annotation import BaseAnnotator
from biocodices.helpers.general import in_groups_of


class DbSNP(BaseAnnotator):
    def annotate(self, rs_list, use_cache=True, use_web=True, parallel=5,
                 sleep_time=10):
        """
        Annotate a list of rs IDs (also accepts a single rs ID).
        When use_cache is True, it will prioritize using Redis cache to get
        info from the given rs. When use_web is True, it will get info from
        dbSNP web. The priority is on the cache, unless explicitely
        inactivated. It returns a dict where the keys are the passed rs IDs.

        * parallel: processes to spawn in parallel when querying dbSNP web.
        * sleep_time: time to sleep between queries.

        Example:
            > dbsnp.annotate(['rs12345', 'rs234'])
            # => {'rs12345': { ... }, 'rs234': { ... }}
        """
        if type(rs_list) == str:
            rs_list = [rs_list]

        rs_list = set(self.remove_non_rs_ids(rs_list))

        info_dict = {}

        if use_cache:
            info_dict.update(self._cache_get(rs_list))
            rs_list = rs_list - info_dict.keys()

        if use_web:
            info_from_api = self._batch_query(rs_list, parallel, sleep_time)
            info_dict.update(info_from_api)

        return info_dict

    def genes(self, rs, use_web=False):
        """Annotate the genes for a given rs."""
        ann = self.annotate(rs, use_web=use_web).get(rs)
        if not ann:
            return []

        mappings = ann.get('assembly', {}).values()
        mappings = chain(*mappings)

        gene_models = [m['geneModel'] for m in mappings if 'geneModel' in m]
        gene_models = chain(*gene_models)

        unique_genes = set()
        for gene in gene_models:
            gene_str = gene['geneSymbol'] or ''
            genes_set = set(gene_str.split('|'))
            unique_genes.update(genes_set)

        return sorted(list(unique_genes))

    @staticmethod
    def _key(rs):
        return 'dbsnp:%s' % rs

    @staticmethod
    def _url(rs):
        path = 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?rs={0}'
        return path.format(rs)

    def _query(self, rs):
        """
        Query NCBI's dbSNP site for a given rs ID. Returns a dict or None.
        """
        url = self._url(rs)
        headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        # print('%s : Query %s' % (rs, url))
        response = requests.get(url, headers)
        # print(' -> %s %s' % (response.status_code, response.reason))

        if response.ok:
            self._cache_set({rs: response.json()})
            return self._cache_get([rs])[rs]
        else:
            return

    def _batch_query(self, rs_list, parallel, sleep_time):
        """
        Query NCBI's dbSNP site for a list of rs IDs. It returns a dict with
        mappings and some data. Runs in <parallel> processes and
        sleeps <sleep> seconds between queries.
        """
        with Pool(parallel) as pool:
            info_dict = {}

            for i, rs_group in enumerate(in_groups_of(parallel, rs_list)):
                if i > 0:
                    print(' Sleep %s seconds' % sleep_time)
                    time.sleep(sleep_time)

                results = {}
                print('Query dbSNP for %s rs IDs' % len(rs_group))
                print(' %s ... %s' % (rs_group[0], rs_group[-1]))
                for rs in set(rs_group):
                    results[rs] = pool.apply_async(self._query, (rs,))

                for rs, result in results.items():
                    try:
                        info = result.get(timeout=20)
                        if info:  # Don't save empty dicts
                            info_dict[rs] = info
                    except TimeoutError:
                        print(rs, 'gave a TimeoutError')

        return info_dict
