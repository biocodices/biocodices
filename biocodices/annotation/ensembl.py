import json
import time
import redis
import requests

from biocodices.helpers.general import in_groups_of


class Ensembl:
    # FIXME: this should check if Redis is present in the system!
    # FIXME: redis config should be read from a YML
    _redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

    def __init__(self, build='GRCh37'):
        self.build = build

    def annotate(self, rs_list, use_cache=True, use_web=True, batch_size=250,
                 sleep_time=15):
        """
        Annotate a list of rs IDs (also accepts a single rs ID).
        When use_cache is True, it will prioritize using Redis cache to get
        info from the given rs. When use_web is True, it will get info from
        Ensembl API. The priority is on the cache, unless it is explicitely
        inactivated. It returns a dict where the keys are the passed rs IDs.

        * batch_size: batch size when hitting the API (max=1000).
        * sleep_time: time to sleep between batch API queries.

        Example:
            > ensembl.annotate(['rs12345', 'rs234'])
            # => {'rs12345': { ... }, 'rs234': { ... }}
        """
        if type(rs_list) == str:
            rs_list = [rs_list]

        rs_list = set([rs for rs in rs_list if rs])
        info_dict = {}

        if use_cache:
            info_dict.update(self._cache_get(rs_list))
            # print('Found %s/%s in Ensembl cache' % (len(info_dict), len(rs_list)))
            rs_list = rs_list - info_dict.keys()

        if use_web:
            info_from_api = self._batch_query(rs_list, batch_size, sleep_time)
            info_dict.update(info_from_api)

        return info_dict


    def _key(self, rs):
        return 'ensembl:%s:%s' % (self.build.lower(), rs)

    def _cache_get(self, rs_list):
        """Get a list of rs IDs Ensembl data from cache."""
        ret = {}
        for rs in rs_list:
            info = self._redis_client.get(self._key(rs))
            if info:
                info = json.loads(info.decode('utf8'))
                ret.update({rs: info})
        return ret

    def _cache_set(self, rs_info_dict, expire_time=None):
        """
        Set the cache for a list of rs IDs. Expects a dict with the form:
        {'rs123': <dict with info about rs123>,
         'rs234': <dict with info about rs234>,
          ... }
        """
        expire_time = expire_time or (60 * 60 * 24 * 30)  # One month

        # Remove empty dicts
        rs_info_dict = {k: v for k, v in rs_info_dict.items() if v}

        # print(' Setting cache for %s keys' % len(rs_info_dict))
        for rs, info in rs_info_dict.items():
            info = json.dumps(info)
            self._redis_client.setex(self._key(rs), expire_time, info)

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
                # print(' Sleep %s seconds' % sleep_time)
                time.sleep(sleep_time)

            # print('Query Ensembl for %s IDs' % len(rs_group))
            # print(' %s ... %s' % (rs_group[0], rs_group[-1]))
            payload = json.dumps({'ids': rs_group, 'phenotypes': '1'})
            response = requests.post(url, headers=headers, data=payload)

            if response.ok:
                ret.update(response.json())
                self._cache_set(response.json())
            else:
                reset_time = int(response.headers['X-RateLimit-Reset'])
                # print(' 400 not OK. Sleeping %s seconds...' % reset_time)
                # print(' Will retry after sleep.')
                time.sleep(reset_time)
                self._batch_query(rs_list)
                # response.raise_for_status()

        return ret
