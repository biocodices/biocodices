import requests
import sys
from myvariant import MyVariantInfo

from biocodices.web_fetchers import MyvariantParser, EnsembleParser


class VariantAnnotator:
    def __init__(self, rs_id):
        self.rs_id = rs_id

    def run(self):
        self.myvariant_series = self._query_myvariant()
        self.ensemble_dict = self._query_ensemble()
        self._parse()

    def _parse(self):
        myvariant_var = MyvariantParser.variant(self.myvariant_series)
        ensemble_var = EnsembleParser.variant(self.ensemble_dict)
        self.variant = {**myvariant_var, **ensemble_var}

        myvariant_pubs = MyvariantParser.publications(self.myvariant_series)
        ensemble_pubs = EnsembleParser.publications(self.ensemble_dict)
        self.publications = myvariant_pubs + ensemble_pubs

    def _query_myvariant(self):
        fields = ['dbsnp', 'dbnsfp', 'grasp', 'gwassnps']
        mv = MyVariantInfo()
        df = mv.query(self.rs_id, fields=fields, as_dataframe=True)
        df.set_index('dbsnp.rsid', inplace=True)

        if len(df.index) > 1:
            msg = 'MyVariant.info returned many entries for {}: {}'
            raise Exception(msg.format(self.rs_id, df.index))

        return df.iloc[0]

    def _query_ensemble(self):
        server = "http://rest.ensembl.org"
        ext = "/variation/human/{}?phenotypes=1".format(self.rs_id)

        headers = { "Content-Type" : "application/json"}
        r = requests.get(server + ext, headers=headers)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        return r.json()
