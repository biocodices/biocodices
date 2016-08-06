import unittest

from biocodices.annotation import Ensembl, DbSNP


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.rs_list = ['rs10192369', 'rs10207628']
        self.dbsnp = DbSNP()
        self.ensembl = Ensembl(build='GRCh37')

    def test_dbsnp_annotate_web(self):
        web_annotation = self.dbsnp.annotate(self.rs_list, use_cache=False,
                                             use_web=True)

        for rs in self.rs_list:
            self.assertIn(rs, web_annotation)

    def test_dbsnp_annotate_cache(self):
        # Make sure annotations are first cached if they don't exist:
        self.dbsnp.annotate(self.rs_list)

        for rs in self.rs_list:
            cache = self.dbsnp.annotate(self.rs_list, use_cache=True,
                                        use_web=False)
            self.assertIn(rs, cache)

    def test_ensembl_annotate_web(self):
        web_annotation = self.ensembl.annotate(self.rs_list, use_cache=False,
                                               use_web=True)

        for rs in self.rs_list:
            self.assertIn(rs, web_annotation)

    def test_ensembl_annotate_cache(self):
        # Make sure annotations are first cached if they don't exist:
        self.ensembl.annotate(self.rs_list)

        for rs in self.rs_list:
            cache = self.ensembl.annotate(self.rs_list, use_cache=True,
                                          use_web=False)
            self.assertIn(rs, cache)

