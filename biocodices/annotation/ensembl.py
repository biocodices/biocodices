import json
import time
import redis
import requests

from biocodices.annotation import AnnotatorWithCache
from biocodices.helpers.general import in_groups_of


class Ensembl(AnnotatorWithCache):
    # FIXME: this should check if Redis is present in the system!
    # FIXME: redis config should be read from a YML
    _redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

    def __init__(self, build='GRCh37'):
        self.build = build

    def genes(self, rs, use_web=False):
        ann = self.annotate(rs, use_web=use_web).get(rs)
        if not ann:
            return []

        gene_names = set()
        for pheno in ann.get('phenotypes', []):
            if 'genes' in pheno and pheno['genes']:
                genes_str = pheno['genes'] or ''
                genes_set = set(genes_str.split(','))
                gene_names.update(genes_set)

        return sorted(list(gene_names))

    def _key(self, rs):
        return 'ensembl:%s:%s' % (self.build.lower(), rs)

    def _batch_query(self, rs_list, batch_size=250, sleep_time=15):
        """
        Takes a list of rs IDs and queries Ensembl via a POST request in batch.
        Returns a dictionary with the rs IDs as keys.
        """
        url = 'http://rest.ensembl.org/variation/homo_sapiens/?phenotypes=1'
        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json'}

        if self.build == 'GRCh37':
            url = url.replace('rest.', 'grch37.rest.')

        # Ensembl API won't take goups of > 1000 identifiers
        if batch_size > 1000:
            batch_size = 1000

        ret = {}
        for i, rs_group in enumerate(in_groups_of(batch_size, rs_list)):
            if i > 0:
                print(' Sleep %s seconds' % sleep_time)
                time.sleep(sleep_time)

            print('Query Ensembl for %s IDs' % len(rs_group))
            print(url)
            print(' %s ... %s' % (rs_group[0], rs_group[-1]))
            payload = json.dumps({'ids': rs_group, 'phenotypes': '1'})
            response = requests.post(url, headers=headers, data=payload)

            if response.ok:
                self._cache_set(response.json())
                ret.update(self._cache_get(rs_group))
            else:
                reset_time = int(response.headers['X-RateLimit-Reset'])
                print(' 400 not OK. Sleeping %s seconds...' % reset_time)
                print(' Will retry after sleep.')
                time.sleep(reset_time)
                self._batch_query(rs_list)
                # response.raise_for_status()

        return ret
