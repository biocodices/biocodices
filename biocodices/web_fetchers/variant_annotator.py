import requests
import sys
from myvariant import MyVariantInfo

from biocodices.web_fetchers import MyvariantParser, EnsembleParser


class VariantAnnotator:
    def run(self, rs):
        """
        Query MyVariant.info and Ensemble for info about an rs ID.
        Returns an Annotation instance that has summary and publications.
        """
        annotation = {}

        annotation['myvariant_data'] = self._query_myvariant(rs)
        annotation['ensemble_data'] = self._query_ensemble(rs)
        variant, publications = self._parse_fetched_data(
            myvariant_df=annotation.myvariant_data,
            ensemble_dict=annotation.ensemble_data,
        )
        annotation['summary'] = variant
        annotation['publications'] = publications

        return annotation

    @staticmethod
    def _parse_fetched_data(myvariant_df, ensemble_dict):
        variant_df = MyvariantParser.parse_annotations(myvariant_df)
        # EnsembleParser.variant(self.ensemble_dict)

        myvariant_pubs = MyvariantParser.publications(myvariant_df)
        ensemble_pubs = EnsembleParser.publications(ensemble_dict)
        publications = myvariant_pubs + ensemble_pubs

        variant_df['publications'] = len(publications)
        for key in ['ancestral_allele', 'most_severe_consequence']:
            variant_df[key] = ensemble_dict[key]

        return variant_df, publications

    def _query_myvariant(self, rs):
        fields = ['all']
        mv = MyVariantInfo()
        df = mv.query(rs, fields=fields, as_dataframe=True)
        df.set_index(['dbsnp.rsid', '_id'], inplace=True)
        return df

    def _query_ensemble(self, rs):
        server = "http://rest.ensembl.org"
        ext = "/variation/human/{}?phenotypes=1".format(rs)

        headers = {"Content-Type": "application/json"}
        r = requests.get(server + ext, headers=headers)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        return r.json()
