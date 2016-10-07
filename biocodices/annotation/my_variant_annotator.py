from multiprocessing import Pool

from myvariant import MyVariantInfo

from biocodices.annotation import BaseAnnotator
from biocodices.helpers.general import in_groups_of


class MyVariantAnnotator(BaseAnnotator):
    def __init__(self):
        self.mv = MyVariantInfo()

    def annotate(self, identifiers, parallel=10, sleep_time=10, use_cache=True,
                 use_web=True):
        """
        Query MyVariant.info for a list of rs IDs or strings like
        'chr1:1234:1235'. Runs in parallel.
        """
        if type(identifiers) == str:
            identifiers = [identifiers]

        identifiers = set(identifiers)
        info_dict = {}

        if use_cache:
            info_dict.update(self._cache_get(identifiers))
            identifiers = identifiers - info_dict.keys()

        if use_web:
            info_from_api = self._batch_query(identifiers, parallel, sleep_time)
            info_dict.update(info_from_api)

        return info_dict

    def _batch_query(self, identifiers, parallel, sleep_time):
        info_dict = {}

        with Pool(parallel) as pool:
            for i, ids_group in enumerate(in_groups_of(parallel, identifiers)):
                if i > 0:
                    print(' Sleep %s seconds' % sleep_time)
                    time.sleep(sleep_time)

                results = {}
                print('Query MyVariant.info for %s IDs' % len(ids_group))
                print(' %s ... %s' % (ids_group[0], ids_group[-1]))
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

    #  def annotate(self, self, identifier, fields=['all']):
        #  """
        #  Query MyVariant.info for info about an rs ID. It also accepts
        #  identifiers like 'chr1:1234-1235'. Returns a dict with the keys:
        #  'myvariant': a pandas DF with selected info from MyVariant.info.
        #  'publications': a list of publications that mention this rs/identifier
                        #  taken both from GRASP via Myvariant.
        #  """
        #  #  myvariant_df, myvariant_publications = \
            #  #  self.parse_myvariant(self.query_myvariant(identifier))

        #  annotation = self.query_myvariant(identifier)

        #  if myvariant_df.empty:
            #  print('No info from MyVariant.info for %s' % identifier)

        #  publications = myvariant_publications

        #  annotation['id'] = identifier

        #  return {
            #  'query': identifier,
            #  'detail': annotation,
            #  'publications': publications
        #  }

    #  @staticmethod
    #  def parse_myvariant(results):
        #  summary = MyvariantParser.parse_query_results(results)
        #  publications = MyvariantParser.parse_publications(results)
        #  return summary, publications

