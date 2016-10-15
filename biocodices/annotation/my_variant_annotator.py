from myvariant import MyVariantInfo

from biocodices.annotation import AnnotatorWithCache


class MyVariantAnnotator(AnnotatorWithCache):
    def __init__(self):
        super(self.__class__, self).__init__()
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

