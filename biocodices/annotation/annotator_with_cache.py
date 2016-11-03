import json
import redis
import time
from multiprocessing import Pool
from biocodices.helpers import in_groups_of

class AnnotatorWithCache():
    """
    Abstract Base Class for Annotators like DbSNP and MyVariant.
    The point of this class is to provide a shared logic of caching the
    responses of remote APIs and of parallelizing requests.
    To use this class, create a new annotator class that implements
    `_query()` and `_key()` methods. Optionally, you can override
    `_batch_query()` for services that you want to parallelize in a different
    way.
    """
    def set_redis_client(self, host, port=6379, db=0):
        klass = self.__class__
        klass._redis_client = redis.StrictRedis(host=host, port=port, db=db)
        conn = klass._redis_client.connection_pool.connection_kwargs
        print('{0} connected to Redis cache: {1}:{2} db={3}'.format(
              self.name, conn['host'], conn['port'], conn['db']))

    def __init__(self, redis_host='localhost', redis_port=6379, redis_db=0):
        self.name = self.__class__.__name__
        # Time for cache to expire:
        self.expire_time = 60 * 60 * 24 * 30 * 5  # Five months in seconds
        self.set_redis_client(redis_host, redis_port, redis_db)

    def annotate(self, ids, parallel=10, sleep_time=10, use_cache=True,
                 use_web=True):
        """
        Annotate a list of variant identifiers, such as rs# IDs or HGVS
        notation. Responses are cached.
        """
        if isinstance(ids, str):
            ids = [ids]
        ids = set(ids)

        info_dict = {}

        if use_cache:
            info_dict.update(self._cache_get(ids))
            ids = ids - info_dict.keys()
        else:
            print(' Not using cache')

        if use_web:
            if ids:
                info_from_api = self._batch_query(ids, parallel, sleep_time)
                info_dict.update(info_from_api)
                ids = ids - info_dict.keys()
        else:
            print(' Not using web')

        if ids:
            print(" No info for %s IDs" % len(ids))

        return info_dict

    def _batch_query(self, ids, parallel, sleep_time):
        print('🌎  Get %s data for %s ids' % (self.name, len(ids)))
        with Pool(parallel) as pool:
            for i, ids_group in enumerate(in_groups_of(parallel, ids)):
                if i > 0:
                    print('  Sleep for %s seconds' % sleep_time)
                    time.sleep(sleep_time)


                print('[{}-{}/{}] {} fetch {} IDs: {} .. {}'.format(
                    (i*parallel)+1, (i*parallel)+parallel, len(ids), self.name,
                    len(ids_group), ids_group[0], ids_group[-1]))

                group_results = {id_: pool.apply_async(self._query, (id_, ))
                                 for id_ in ids_group}

                for id_, result in group_results.items():
                    group_results[id_] = result.get()

                print('  Set cache for %s IDs' % len(ids_group))
                self._cache_set(group_results)

        # After caching all the responses, use the logic in #_cache_get()
        # to bring them from cache. The jsons will be json-loaded, the xml
        # will be left as they are, etc.
        return self._cache_get(ids)

    def _cache_set(self, info_dict):
        """
        Set the cache for a list of IDs. Expects a dict with the form:
        {
            identifier1: <dict with info about the id1> or data_string,
            identifier2: <dict with info about id2> or data_string,
            ...
        }
        It will json-dump the dicts.
        """
        # Remove empty responses
        info_dict = {k: v for k, v in info_dict.items() if v}

        for identifier, info in info_dict.items():
            if isinstance(info, dict):
                info = json.dumps(info)
            self._redis_client.setex(self._key(identifier), self.expire_time, info)

    def _cache_get(self, ids):
        """
        Get a list of IDs data from cache. Returns a dict with the form:
        {
            identifier1: <dict with info about the id1> or data_string,
            identifier2: <dict with info about id2> or data_string,
            ...
        }
        """
        info_dict = {}
        for identifier in ids:
            cached_info = self._redis_client.get(self._key(identifier))
            if cached_info:
                try:
                    cached_info = json.loads(cached_info.decode('utf8'))
                except ValueError:  # Not valid JSON, leave response as is
                    cached_info = cached_info.decode('utf8')

                info_dict[identifier] = cached_info

        if len(ids) > 1:
            msg = '📂 %s found %s/%s in cache'
            print(msg % (self.name, len(info_dict), len(ids)))

        return info_dict

    def _query(self):
        raise NotImplementedError()

    def _key(self):
        raise NotImplementedError()

