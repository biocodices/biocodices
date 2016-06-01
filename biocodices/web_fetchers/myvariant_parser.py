class MyvariantParser:
    @classmethod
    def variant(cls, myvariant_series):
        var_dict = {}
        functions = [cls._parse_dbsnp_fields,
                     cls._parse_dbnsfp_fields,
                     cls._parse_grasp_fields,
                     cls._parse_gwassnps_fields]

        for function in functions:
            var_dict.update(function(myvariant_series))

        return var_dict

    @classmethod
    def publications(cls, myvariant_series):
        publications = []
        publications += cls._parse_grasp_publications_field(myvariant_series)
        return publications

    @staticmethod
    def _parse_grasp_publications_field(s):
        if len(s.filter(regex='grasp').index) == 0:
            return []

        publications = []
        for pub in s['grasp.publication']:
            publications.append({
                'pubmed_id': pub['pmid'],
                'title': pub['title'],
                'phenotype': ''.join(pub['phenotype']),
                # ^ applies to str and lists
                'date': pub['date_pub'],
            })

        return publications

    @staticmethod
    def _parse_dbsnp_fields(s):
        if len(s.filter(regex='dbsnp').index) == 0:
            return {}

        return {
            'name': s.name,
            'source': 'dbsnp',
            'hg19_start': s['dbsnp.hg19.start'],
            'hg19_end': s['dbsnp.hg19.end'],
            'chromosome': s['dbsnp.chrom'],
            'alleles': [d['allele'] for d in s['dbsnp.alleles']],
            'ref': s['dbsnp.ref'],
            'alt': s['dbsnp.alt'],
        }

    @staticmethod
    def _parse_grasp_fields(s):
        if len(s.filter(regex='grasp').index) == 0:
            return {}

        return {
            'polyphen2': s.get('grasp.polyphen2'),
            'sift': s.get('grasp.sift'),
            'mapped_gene_grasp': s.get('grasp.in_gene'),
        }

    @staticmethod
    def _parse_gwassnps_fields(s):
        if len(s.filter(regex='gwassnps').index) == 0:
            return {}

        return {
            'gwassnps_pubmed': s['gwassnps.pubmed'],
            'gwassnps_pvalue': s['gwassnps.pvalue'],
            'risk_allele': s['gwassnps.risk_allele'],
            'gwassnps_title': s['gwassnps.title'],
            'gwassnps_trait': s['gwassnps.trait'],
        }

    @staticmethod
    def _parse_dbnsfp_fields(s):
        if len(s.filter(regex='dbnsfp').index) == 0:
            return {}

        return {
            'ancestral_allele': s['dbnsfp.ancestral_allele'].upper(),
            'strand': s['dbnsfp.cds_strand'],
            'polyphen2_hdiv_prediction': s['dbnsfp.polyphen2.hdiv.pred'],
            'polyphen2_hdiv_score': s['dbnsfp.polyphen2.hdiv.score'],
            'polyphen2_hdiv_rankscore': s['dbnsfp.polyphen2.hdiv.rankscore'],
            'polyphen2_hvar_prediction': s['dbnsfp.polyphen2.hvar.pred'],
            'polyphen2_hvar_score': s['dbnsfp.polyphen2.hvar.score'],
            'polyphen2_hvar_rankscore': s['dbnsfp.polyphen2.hvar.rankscore'],
            'sift_prediction': s['dbnsfp.sift.pred'],
        }
