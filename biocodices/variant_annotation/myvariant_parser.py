import pandas as pd
import re


class MyvariantParser:
    @classmethod
    def parse_query_results(cls, results):
        """Deals with MyVariant.info restuls dictionary. Returns a pandas
        DataFrame with some selected fields and one entry per HGVS id."""
        entries = []

        if results['total'] > 0:
            for hit in results['hits']:
                entry = {'hgvs_id': hit['_id']}
                entry.update(cls.parse_dbsnp(hit.get('dbsnp')))
                entry.update(cls.parse_cadd(hit.get('cadd')))
                entry.update(cls.parse_wellderly(hit.get('wellderly')))
                entry.update(cls.parse_snpeff(hit.get('snpeff')))
                entry.update(cls.parse_grasp(hit.get('grasp')))
                entries.append(entry)

        df = pd.DataFrame(entries)

        # Generate a list of unique genes mentioned across fields
        uniq_genes = df.filter(regex='gene').sum(axis=1).map(lambda l: set(l))
        df['genes'] = uniq_genes.map(list)

        return df

    @classmethod
    def parse_publications(cls, results):
        publications = []

        if results['total'] > 0:
            for hit in results['hits']:
                publications += cls.grasp_publications(hit.get('grasp'))

        return publications

    @classmethod
    def parse_grasp(cls, grasp_dict):
        """Return a flat dictionary with some GRASP fields."""
        grasp_dict = grasp_dict or {}
        summary = {}

        summary['grasp_ancestries'] = \
            grasp_dict.get('gwas_ancestry_description')

        summary['grasp_genes'] = []
        summary['grasp_pubmed_ids'] = [pub['pubmed_id'] for pub in
                                       cls.grasp_publications(grasp_dict)
                                       if pub['pubmed_id']]

        gene_names_str = grasp_dict.get('in_gene')
        if gene_names_str:
            # Example entries: '(FAM47E)(FAM47E-STBD1)', '(PARK7)'
            gene_names = re.findall('\((.+?)\)', gene_names_str)
            summary['grasp_genes'] += gene_names

        return summary

    @staticmethod
    def listify(element):
        """Transform a datum that might not be a list into a list of one
        element, but leave the item as it is if it's already a list.
        Return an empty list if the passed element is None."""
        if element is None:
            return []
        if type(element) == list:
            return element
        return [element]

    @classmethod
    def grasp_publications(cls, grasp_dict):
        publications_summary = []

        if not grasp_dict or 'publication' not in grasp_dict:
            return publications_summary

        for pub in cls.listify(grasp_dict['publication']):
            publications_summary.append({
                'pubmed_id': pub.get('pmid'),
                'title': pub.get('title'),
                'phenotype': ''.join(pub.get('phenotype')),
                # ^ The 'phenotype' field of publications is sometimes a
                # sentence split in two forming a list (!).
                # This join operates on that, fixing it, and it ALSO operates
                # on regular string fields, in that case leaving the string
                # as it is.
                'date': pub.get('date_pub'),
            })

        return publications_summary

    @staticmethod
    def parse_dbsnp(dbsnp_dict):
        """Return a flat dictionary with some dbSNP fields."""
        dbsnp_dict = dbsnp_dict or {}
        summary = {}

        simple_fields = ['alt', 'ref', 'chrom', 'class']
        for field in simple_fields:
            summary['dbsnp_' + field] = dbsnp_dict.get(field)

        hg19 = dbsnp_dict.get('hg19')
        if hg19:
            summary['dbsnp_hg19_start'] = hg19.get('start')
            summary['dbsnp_hg19_end'] = hg19.get('end')

        return summary

    @classmethod
    def parse_snpeff(cls, snpeff_dict):
        """Return a flat dictionary with some SnpEff fields."""
        snpeff_dict = snpeff_dict or {}
        summary = {}

        unique_genes = set()
        summary['snpeff_ann'] = []

        for annotation in cls.listify(snpeff_dict.get('ann')):
            gene_name = annotation.get('gene_name')
            unique_genes.add(gene_name)

            summary_annotation = (
                annotation.get('effect'),
                gene_name,
                annotation.get('putative_impact'),
            )  # Tuple with some selected fields
            summary['snpeff_ann'].append(summary_annotation)

        summary['snpeff_genes'] = list(unique_genes)

        return summary

    @classmethod
    def parse_wellderly(cls, wellderly_dict):
        """Return a flat dictionary with some Wellderly fields."""
        wellderly_dict = wellderly_dict or {}
        summary = {}

        summary['wellderly_genes'] = cls.listify(wellderly_dict.get('gene'))

        return summary

    @classmethod
    def parse_cadd(cls, cadd_dict):
        """Return a flat dictionary with some CADD fields."""
        cadd_dict = cadd_dict or {}
        summary = {}

        simple_fields = ['annotype', 'consequence']
        for field in simple_fields:
            summary['cadd_' + field] = cadd_dict.get(field)

        summary['cadd_genes'] = []
        for feature in cls.listify(cadd_dict.get('gene')):
            # Example result:
            #
            #  [{'feature_id': 'ENST00000451424',
            #    'gene_id': 'ENSG00000117242',
            #    'genename': 'PINK1-AS'},
            #   {'ccds_id': 'CCDS211.1',
            #    'feature_id': 'ENST00000321556',
            #    'gene_id': 'ENSG00000158828',
            #    'genename': 'PINK1'}]
            gene_name = feature.get('genename')
            if gene_name:
                summary['cadd_genes'].append(gene_name)

        return summary
