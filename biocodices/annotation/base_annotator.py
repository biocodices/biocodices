import json
import redis


class BaseAnnotator():
    _redis_client = redis.StrictRedis(host='localhost', port=6379, db=0)

    def _cache_set(self, id_info_dict, expire_time=None):
        """
        Set the cache for a list of IDs. Expects a dict with the form:
        {'rs123': <dict with info about rs123>,
         'rs234': <dict with info about rs234>,
          ... }
        It will json-dump the dicts.
        """
        expire_time = expire_time or (60 * 60 * 24 * 30 * 5)  # Five months

        # Remove empty dicts
        id_info_dict = {k: v for k, v in id_info_dict.items() if v}

        # print(' Setting cache for %s keys' % len(id_info_dict))
        for identifier, info in id_info_dict.items():
            info = json.dumps(info)
            self._redis_client.setex(self._key(identifier), expire_time, info)

    def _cache_get(self, id_list):
        """Get a list of rs IDs Ensembl data from cache."""
        ret = {}
        for identifier in id_list:
            info = self._redis_client.get(self._key(identifier))
            if info:
                info = json.loads(info.decode('utf8'))
                ret.update({identifier: info})

        if len(id_list) > 1:
            print('Found %s/%s in dbSNP cache' % (len(ret), len(id_list)))
        return ret

    def annotate(self):
        raise NotImplementedError()

    def _query(self):
        raise NotImplementedError()

    def _key(self):
        raise NotImplementedError()

    @staticmethod
    def remove_non_rs_ids(id_list):
        non_rs_ids = [id_ for id_ in id_list if id_ and 'rs' not in id_]

        if non_rs_ids:
            print('Leave out %s ids without "rs"' % len(non_rs_ids))

        return [rs for rs in id_list if rs and rs not in non_rs_ids]
