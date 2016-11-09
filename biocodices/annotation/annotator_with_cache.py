import sys
import ujson
import redis
import time
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm
from biocodices.helpers import in_groups_of

class AnnotatorWithCache():
    """
    Abstract Base Class for Annotators like DbSNP and MyVariant.
    The point of this class is to provide a shared logic of caching the
    responses of remote APIs and of parallelizing requests.
    To use this class, create a new annotator class that implements
    `_query()` and `_key()` methods. Optionally, you can override
    `_batch_query()` for services that you want to parallelize in a different
    way or not paralleliza at all.
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
        self.set_redis_client(redis_host, redis_port, redis_db)

    def annotate(self, ids, parallel=10, sleep_time=10, use_cache=True,
                 use_web=True):
        """
        Annotate a list of variant identifiers, such as rs# IDs or HGVS
        notation. Responses are cached.
        """
        if isinstance(ids, str) or isinstance(ids, int):
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

    def _query_and_set_cache(self, id_):
        response = self._query(id_)
        if response:
            self._cache_set({id_: response})
        return response

    def _batch_query(self, ids, parallel, sleep_time):
        grouped_ids = list(in_groups_of(parallel, ids))

        msg = '🌎 Get {} {} entries in {} batches ({} items/batch)'
        msg += ', sleeping {}s between batches'
        print(msg.format(len(ids), self.name, len(grouped_ids), parallel, sleep_time))

        with ThreadPoolExecutor(max_workers=parallel) as executor:
            sys.stdout.flush()  # Necessary for the progress bar correct display
            for i, ids_group in enumerate(tqdm(grouped_ids, total=len(grouped_ids))):
                if i > 0:
                    time.sleep(sleep_time)

                executor.map(self._query_and_set_cache, ids_group)

        # After caching all the responses, use the logic in #_cache_get()
        # to bring them from cache. The jsons will be json-loaded, the xml
        # will be left as they are, etc.
        return self._cache_get(ids, verbose=False)

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
        data_to_cache = {}
        for id_, value in info_dict.items():
            if not value:
                continue
            if isinstance(value, dict) or isinstance(value, list):
                value = ujson.dumps(value)
            data_to_cache[self._key(id_)] = value

        self._redis_client.mset(data_to_cache)

    def _cache_get(self, ids, verbose=True):
        """
        Get a list of IDs data from cache. Returns a dict with the form:
        {
            identifier1: <dict with info about the id1> or data_string,
            identifier2: <dict with info about id2> or data_string,
            ...
        }
        """
        keys = [self._key(id_) for id_ in ids]
        info_dict = {k: v for k, v in zip(ids, self._redis_client.mget(keys)) if v}

        for k, v in info_dict.items():
            v = self._decode_and_try_deserialize(v)
            info_dict[k] = v

        if verbose and len(ids) > 1:
            msg = '📂 %s found %s/%s in cache'
            print(msg % (self.name, len(info_dict), len(ids)))

        return info_dict

    def _query(self):
        raise NotImplementedError()

    def _key(self):
        raise NotImplementedError()

    @staticmethod
    def _decode_and_try_deserialize(data):
        try:
            return ujson.loads(data.decode('utf8'))
        except ValueError:  # Not valid JSON, leave response as is
            return data.decode('utf8')

