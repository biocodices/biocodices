import json
import time
import redis
import requests

from biocodices.helpers.general import in_groups_of


class Ensembl:
    # FIXME: this should check if Redis is present in the system!
    # FIXME: redis config should be read from a YML
    redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

    @staticmethod
    def key(rs):
        return 'ensembl:%s' % rs

    @classmethod
    def cache_get(cls, rs_list):
        """Get a list of rs IDs Ensembl data from cache."""
        ret = {}
        for rs in rs_list:
            info = cls.redis_client.get(cls.key(rs))
            if info:
                info = json.loads(info.decode('utf8'))
                ret.update({rs: info})
        return ret

    @classmethod
    def cache_set(cls, rs_info_dict, expire_time=None):
        """
        Set the cache for a list of rs IDs. Expects a dict with the form:
        {'rs123': <dict with info about rs123>,
         'rs234': <dict with info about rs234>,
          ... }
        """
        expire_time = expire_time or (60 * 60 * 24 * 30)  # One month
        rs_info_dict = {k: v for k, v in rs_info_dict.items() if v}
        print('Setting cache for %s keys' % len(rs_info_dict))
        for rs, info in rs_info_dict.items():
            if not info:
                continue  # Don't save empty dicts to cache
            info = json.dumps(info)
            cls.redis_client.setex(cls.key(rs), expire_time, info)

    @classmethod
    def batch_query(cls, rs_list, build='GRCh37', use_cache=True, sleep=3):
        """
        Takes a list of rs IDs and queries Ensembl. Returns a dictionary with
        the rs IDs as keys. It can annotate with builds GRCh37 and GRCh38.
        """
        rs_list = set([rs for rs in rs_list if rs])
        url = 'http://rest.ensembl.org/variation/homo_sapiens?phenotypes=1'
        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json'}

        if build == 'GRCh37':
            url = url.replace('rest.', 'grch37.rest.')

        ret = {}
        if use_cache:
            ret.update(cls.cache_get(rs_list))
            print('Found %s/%s in cache' % (len(ret), len(rs_list)))
            rs_list = rs_list - ret.keys()
            # ^ Remove the annotated ones from cache

        # Ensembl API won't take goups of > 1000 identifiers
        for rs_group in in_groups_of(100, rs_list):
            print('Query Ensembl for %s identifiers...' % len(rs_group))
            payload = json.dumps({'ids': rs_group, 'phenotypes': '1'})
            response = requests.post(url, headers=headers, data=payload)

            if response.ok:
                ret.update(response.json())
                cls.cache_set(response.json())
                print(' -> 200 OK!')
            else:
                reset_time = int(response.headers['X-RateLimit-Reset'])
                print(' -> 400 not OK. Sleeping %s seconds...' % reset_time)
                print('    (Will retry after sleep).')
                time.sleep(reset_time)
                cls.batch_query(rs_list)
                # response.raise_for_status()

            print(' -> Sleep %s seconds' % sleep)
            time.sleep(sleep)

        return ret
