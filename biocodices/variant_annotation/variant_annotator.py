import requests
from multiprocessing import Pool
import pandas as pd
from myvariant import MyVariantInfo

from biocodices.variant_annotation import MyvariantParser, EnsembleParser
from biocodices.helpers.general import restful_api_query


class VariantAnnotator:
    def __init__(self):
        self.full_annotations = {}
        self.annotations = pd.DataFrame({})

    @staticmethod
    def summary_annotation(rs):
        return VariantAnnotator().annotate(rs)['summary']

    def annotate_in_batch(self, rs_list, processes=10, timeout=10):
        """Query MyVariant.info and Ensemble for a list of rs IDs. Runs in
        parallel. Returns a merged pandas DataFrame with the SNPs as
        indices."""
        with Pool(processes) as pool:
            results = [pool.apply_async(self.summary_annotation, (rs, ))
                       for rs in rs_list]
            annotations = [result.get(timeout=timeout)
                           for result in results]

        return pd.concat(annotations)

    @classmethod
    def annotate(cls, rs):
        """
        Query MyVariant.info and Ensemble for info about an rs ID.
        Returns a dict with the following data:
            'myvariant': a dict with selected info from MyVariant.info
            'ensemble': a dict with the json response from Ensemble
            'summary': a handy summary of the info from the above sources
            'publications': a list of publications that mention this rs
        """
        annotation = {}

        annotation['myvariant'] = cls.query_myvariant(rs)

        #  annotation['ensemble'] = cls.query_ensemble(rs)
        #  variant, publications = self._parse_fetched_data(
            #  myvariant_df=annotation['myvariant'],
            #  ensemble_dict=annotation['ensemble'],
        #  )
        #  annotation['summary'] = variant
        #  annotation['publications'] = publications

        #  self.full_annotations[rs] = annotation
        #  self.annotations = self.annotations.append(annotation['summary'])

        return annotation

    @staticmethod
    def _parse_fetched_data(myvariant_df, ensemble_dict):
        variant_df = MyvariantParser.parse_annotations(myvariant_df)

        if ensemble_dict:
            for key in ['ancestral_allele', 'most_severe_consequence']:
                variant_df[key] = ensemble_dict[key]

        myvariant_pubs = MyvariantParser.publications(myvariant_df)
        ensemble_pubs = EnsembleParser.publications(ensemble_dict)
        publications = myvariant_pubs + ensemble_pubs

        variant_df['publications'] = len(publications)

        return variant_df, publications

    @staticmethod
    def query_myvariant(rs):
        results = MyVariantInfo().query(rs, fields=['all'])
        df = MyvariantParser.parse_query_results(results)
        df['rs_id'] = rs
        return df

    @staticmethod
    def query_ensemble(rs):
        url = "http://rest.ensembl.org"
        url += "/variation/human/{}?phenotypes=1".format(rs)
        try:
            return restful_api_query(url)
        except requests.exceptions.HTTPError as error:
            print(error)
            return {}

    @staticmethod
    def query_cellbase(chromosome, position, allele, key='snp_phenotype'):
        """Queries CellBase restful API. Default key 'snp_phenotype' will fetch
        info about the pheno effect of a SNP. Other options are
        'mutation_phenotype' for COSMIC database and 'consequence_type'."""
        url = 'http://ws.bioinfo.cipf.es/cellbase/rest/latest/hsa/genomic/variant/'
        url += '{chromosome}:{position}:{allele}'
        url += '/{key}?of=json'

        return restful_api_query(url)
