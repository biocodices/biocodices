class EnsembleParser:
    @staticmethod
    def variant(ensemble_dict):
        variant = {
            'MAF': ensemble_dict['MAF'],
            'ensemble_ancestral_allele': ensemble_dict['ancestral_allele'],
            'ensemble_minor_allele': ensemble_dict['minor_allele'],
            'ensemble_most_severe_consequence': ensemble_dict['most_severe_consequence'],
        }

        for mapping in ensemble_dict['mappings']:
            if mapping['assembly_name'] == 'GRCh38':
                    variant['GRCh38_start'] = mapping['start']
                    variant['GRCh38_end'] = mapping['end']

        return variant

    @staticmethod
    def publications(ensemble_dict):
        publications = []

        for study in ensemble_dict['phenotypes']:
            if study['source'] == 'ClinVar':
                # ClinVar via Ensemble doesn't include publication info
                continue

            print(study)
            publication = {
                'genes': study['genes'],
                'pvalue': study['pvalue'],
                'trait': study['trait'],
            }
            if 'study' in study:
                source, pub_id = study['study'].split(':')
                if source == 'PMID':
                    publication['pubmed_id'] = pub_id
                if source == 'MIM':
                    publication['omim_id'] = pub_id

            publications.append(publication)

        return publications
