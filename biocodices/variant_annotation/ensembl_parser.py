import pandas as pd


class EnsemblParser:
    @staticmethod
    def parse_query_results(ensembl_dict):
        summary = {}
        if not ensembl_dict:
            return summary

        simple_fields = ['MAF', 'ancestral_allele', 'minor_allele',
                         'most_severe_consequence']

        for field in simple_fields:
            summary['ensembl_' + field] = ensembl_dict.get(field)

        mappings = ensembl_dict.get('mappings') or []
        for mapping in mappings:
            assembly = mapping['assembly_name']
            summary['{}_start'.format(assembly)] = mapping.get('start')
            summary['{}_end'.format(assembly)] = mapping.get('end')

        return pd.Series(summary).to_frame().transpose()

    @staticmethod
    def parse_publications(ensembl_dict):
        publications = []

        if not ensembl_dict:
            return publications

        phenotypes = ensembl_dict['phenotypes'] or []
        for study in phenotypes:
            if study.get('source') == 'ClinVar':
                # ClinVar via Ensembl doesn't include publication info
                continue

            publication = {k: v for k, v in study.items()
                           if study.get(k) and
                           k in ['genes', 'pvalue', 'trait']}

            study_id = study.get('study')
            if study_id:
                source, pub_id = study_id.split(':')
                if source == 'PMID':
                    publication['pubmed_id'] = pub_id
                if source == 'MIM':
                    publication['omim_id'] = pub_id

            publications.append(publication)

        return publications
