import requests
from multiprocessing import Pool
import pandas as pd
from myvariant import MyVariantInfo

from biocodices.variant_annotation import MyvariantParser, EnsemblParser
from biocodices.helpers.general import restful_api_query


class VariantAnnotator:
    def __init__(self):
        self.full_annotations = {}
        self.annotations = pd.DataFrame({})

    def annotate_in_batch(self, rs_list, processes=10, timeout=10):
        """Query MyVariant.info and Ensembl for a list of rs IDs. Runs in
        parallel. Returns a merged pandas DataFrame with the SNPs as
        indices."""
        with Pool(processes) as pool:
            results = [pool.apply_async(self.annotate, (rs, ))
                       for rs in rs_list]
            results = [result.get(timeout=timeout) for result in results]

        annotations = pd.concat([ann['annotation'] for ann in results],
                                ignore_index=True)
        index = ['dbsnp_chrom', 'dbsnp_hg19_start', 'rs_id', 'hgvs_id']
        annotations.set_index(index, inplace=True)
        publications = {ann['rs']: ann['publications'] for ann in results}
        return {'annotations': annotations, 'publications': publications}

    @classmethod
    def annotate(cls, rs):
        """
        Query MyVariant.info and Ensembl for info about an rs ID.
        Returns a dict with the following keys:
        'myvariant': a pandas DF with selected info from MyVariant.info.
        'ensembl': a dict with selected info from Ensembl.
        'publications': a list of publications that mention this rs
                        taken both from GRASP via Myvariant and Ensembl.
        """
        myvariant_df, myvariant_publications = cls.query_myvariant(rs)
        ensembl_df, ensembl_publications = cls.query_ensembl(rs)

        annotation = myvariant_df.join(ensembl_df)
        annotation['rs_id'] = rs
        publications = myvariant_publications + ensembl_publications
        return {
            'rs': rs,
            'annotation': annotation,
            'publications': publications
        }

    @staticmethod
    def query_myvariant(rs):
        results = MyVariantInfo().query(rs, fields=['all'])
        summary = MyvariantParser.parse_query_results(results)
        publications = MyvariantParser.parse_publications(results)
        return summary, publications

    @staticmethod
    def query_ensembl(rs):
        url = "http://rest.ensembl.org"
        url += "/variation/human/{}?phenotypes=1".format(rs)
        try:
            results = restful_api_query(url)
        except requests.exceptions.HTTPError as error:
            # Ensembl API returns an error 400 when the query gives no results.
            print(error)
            results = {}

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
