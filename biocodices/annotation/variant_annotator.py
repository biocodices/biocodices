import requests
from multiprocessing import Pool
import pandas as pd
from myvariant import MyVariantInfo

from biocodices.annotation import MyvariantParser, EnsemblParser
from biocodices.helpers.general import restful_api_query


class VariantAnnotator:
    @classmethod
    def annotate_in_batch(cls, identifiers, processes=10, timeout=10,
                          myvariant=True, ensembl=True):
        """Dead code, left here not to break old notebooks."""
        with Pool(processes) as pool:
            annotate_args = dict(myvariant=myvariant, ensembl=ensembl)
            results = [pool.apply_async(cls.annotate, (identifier,), annotate_args)
                       for identifier in set(identifiers)]
            results = [result.get(timeout=timeout) for result in results]

        annotations = pd.concat([ann['detail'] for ann in results],
                                ignore_index=True)
        #  index = ['dbsnp_chrom', 'dbsnp_hg19_start', 'rs_id', 'hgvs_id']
        #  annotations.set_index(index, inplace=True)
        publications = {ann['identifier']: ann['publications'] for ann in results}
        return {'detail': annotations, 'publications': publications}

    @classmethod
    def annotate(cls, identifier, myvariant=True, ensembl=True,
                 myvariant_fields=['all']):
        """Dead code, left here not to break old notebooks."""
        if myvariant:
            myvariant_df, myvariant_publications = \
                cls.parse_myvariant(cls.query_myvariant(identifier))
            if myvariant_df.empty:
                print('No info from MyVariant.info for %s' % identifier)
        if ensembl:
            ensembl_df, ensembl_publications = \
                cls.parse_ensembl(cls.query_ensembl(identifier))

        if myvariant and ensembl:
            annotation = myvariant_df.join(ensembl_df)
            publications = myvariant_publications + ensembl_publications
        elif myvariant:
            annotation = myvariant_df
            publications = myvariant_publications
        elif ensembl:
            annotation = ensembl_df
            publications = ensembl_publications

        annotation['rs_id'] = identifier  # Leaving this for backwards compatibility
        annotation['query'] = identifier
        return {
            'identifier': identifier,
            'detail': annotation,
            'publications': publications
        }

    @staticmethod
    def query_myvariant(identifier, fields=['all']):
        """Dead code, left here not to break old notebooks."""
        return MyVariantInfo().query(identifier, fields=fields)

    @staticmethod
    def parse_myvariant(results):
        """Dead code, left here not to break old notebooks."""
        summary = MyvariantParser.parse_query_results(results)
        publications = MyvariantParser.parse_publications(results)
        return summary, publications

    @staticmethod
    def query_ensembl(rs):
        url = "http://grch37.rest.ensembl.org"
        url += "/variation/human/{}?phenotypes=1".format(rs)
        try:
            print('Query', url)
            results = restful_api_query(url)
        except requests.exceptions.HTTPError as error:
            # Ensembl API returns an error 400 when the query gives no results.
            print(error)
            results = {}
        return results

    @staticmethod
    def parse_ensembl(results):
        summary = EnsemblParser.parse_query_results(results)
        publications = EnsemblParser.parse_publications(results)
        return summary, publications

    @staticmethod
    def query_cellbase(chromosome, position, allele, key='snp_phenotype'):
        """Queries CellBase restful API. Default key 'snp_phenotype' will fetch
        info about the pheno effect of a SNP. Other options are
        'mutation_phenotype' for COSMIC database and 'consequence_type'."""
        url = 'http://ws.bioinfo.cipf.es/cellbase/rest/latest/hsa/genomic/variant/'
        url += '{chromosome}:{position}:{allele}'
        url += '/{key}?of=json'

        return restful_api_query(url)
