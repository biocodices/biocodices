from myvariant import MyVariantInfo

from biocodices.web_fetchers import MyvariantParser, EnsembleParser
from biocodices.helpers.general import restful_api_query


class VariantAnnotator:
    def run(self, rs):
        """
        Query MyVariant.info and Ensemble for info about an rs ID.
        Returns an Annotation instance that has summary and publications.
        """
        annotation = {}

        annotation['myvariant_data'] = self.query_myvariant(rs)
        annotation['ensemble_data'] = self.query_ensemble(rs)
        variant, publications = self._parse_fetched_data(
            myvariant_df=annotation['myvariant_data'],
            ensemble_dict=annotation['ensemble_data'],
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

    @staticmethod
    def query_myvariant(rs):
        fields = ['all']
        mv = MyVariantInfo()
        df = mv.query(rs, fields=fields, as_dataframe=True)
        df.set_index(['dbsnp.rsid', '_id'], inplace=True)
        return df

    @staticmethod
    def query_ensemble(rs):
        server = "http://rest.ensembl.org"
        ext = "/variation/human/{}?phenotypes=1".format(rs)
        return restful_api_query(server + ext)

    @staticmethod
    def query_cellbase(chromosome, position, allele, key='snp_phenotype'):
        """Queries CellBase restful API. Default key 'snp_phenotype' will fetch
        info about the pheno effect of a SNP. Other options are
        'mutation_phenotype' for COSMIC database and 'consequence_type'."""
        url = 'http://ws.bioinfo.cipf.es/cellbase/rest/latest/hsa/genomic/variant/'
        url += '{chromosome}:{position}:{allele}'
        url += '/{key}?of=json'

        print(url)
        return restful_api_query(url)
