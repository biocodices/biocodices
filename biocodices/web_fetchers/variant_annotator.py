import requests
import sys
from collections import namedtuple
from myvariant import MyVariantInfo

from biocodices.web_fetchers import MyvariantParser, EnsembleParser


class VariantAnnotator:
    def run(self, rs):
        """Query MyVariant.info and Ensemble for info about an rs ID."""
        return self._parse_fetched_data(
            myvariant_df=self._query_myvariant(rs),
            ensemble_dict=self._query_ensemble(rs),
        )

    @staticmethod
    def _parse_fetched_data(myvariant_df, ensemble_dict):
        Info = namedtuple('VariantAnnotation',
                          ['myvariant', 'ensemble', 'publications'])

        variant_df = MyvariantParser.parse_annotations(myvariant_df)
        # EnsembleParser.variant(self.ensemble_dict)

        myvariant_pubs = MyvariantParser.publications(myvariant_df)
        ensemble_pubs = EnsembleParser.publications(ensemble_dict)
        publications = myvariant_pubs + ensemble_pubs

        return Info(variant_df, ensemble_dict, publications)

    def _query_myvariant(self, rs):
        # fields = ['dbsnp', 'dbnsfp', 'grasp', 'gwassnps']
        fields = ['all']
        mv = MyVariantInfo()
        df = mv.query(rs, fields=fields, as_dataframe=True)
        df.set_index(['dbsnp.rsid', '_id'], inplace=True)
        return df

    def _query_ensemble(self, rs):
        server = "http://rest.ensembl.org"
        ext = "/variation/human/{}?phenotypes=1".format(rs)

        headers = { "Content-Type" : "application/json"}
        r = requests.get(server + ext, headers=headers)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        return r.json()
