import time
from multiprocessing import Pool

from myvariant import MyVariantInfo

from biocodices.annotation import AnnotatorWithCache
from biocodices.helpers.general import in_groups_of


class MyVariantAnnotator(AnnotatorWithCache):
    def __init__(self):
        self.mv = MyVariantInfo()

    @staticmethod
    def _key(identifier):
        return 'myvariant:%s' % identifier

    def _query(self, identifier, fields=['all']):
        response = self.mv.query(identifier, fields=fields)
        if response:
            self._cache_set({identifier: response})
            return self._cache_get([identifier])[identifier]
        else:
            return {}

    def _batch_query(self, identifiers, parallel, sleep_time):
        info_dict = {}

        with Pool(parallel) as pool:
            for i, ids_group in enumerate(in_groups_of(parallel, identifiers)):
                if i > 0:
                    # print(' Sleep %s seconds' % sleep_time)
                    time.sleep(sleep_time)

                msg = 'Batch {0}/{1}. Query MyVariant.info for {2} variants'
                print(msg.format(i+1, len(identifiers)//parallel, len(ids_group)))
                print(' %s ... %s' % (ids_group[0], ids_group[-1]))

                results = {}

                for id_ in set(ids_group):
                    results[id_] = pool.apply_async(self._query, (id_,))

                for id_, result in results.items():
                    try:
                        info = result.get(timeout=20)
                        if info:  # Don't save empty dicts
                            info_dict[id_] = info
                    except TimeoutError:
                        print(id_, 'gave a TimeoutError')

        return info_dict

