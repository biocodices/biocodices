import pandas as pd


class MyvariantParser:
    @classmethod
    def publications(cls, myvariant_series):
        publications = []
        publications += cls._parse_grasp_publications_field(myvariant_series)
        return publications

    @classmethod
    def parse_annotations(cls, variant_df):
        gathered_info = {}

        for colname in ['wellderly.adviser_score',
                        'wellderely.gene']:
            if colname in variant_df:
                gathered_info[colname] = variant_df[colname].dropna()

        for colname in ['cadd.gene.genename',
                        'cadd.annotype',
                        'cadd.consequence']:
            if colname in variant_df:
                gathered_info[colname] = variant_df[colname].dropna() \
                                                            .map(cls.listify)

        if 'cadd.gene' in variant_df:
            gathered_info['cadd.gene'] = \
                variant_df['cadd.gene'].dropna().map(cls.cadd_genenames)

        if cls.info_present_from('dbsnp', variant_df):
            gathered_info.update({
                'hg19_start': variant_df['dbsnp.hg19.start'],
                'chrom': variant_df['dbsnp.chrom'],
                'ref': variant_df['dbsnp.ref'],
                'alt': variant_df['dbsnp.alt'],
            })

        if 'snpeff.ann' in variant_df:
            snpeff_annot = variant_df['snpeff.ann'].dropna().map(cls.snpeff_ann)
            gathered_info['snpeff_annotation'] = snpeff_annot

        annotations = pd.DataFrame(gathered_info)

        if cls.info_present_from('dbsnp', variant_df):
            var_allele = annotations.index.get_level_values('_id') \
                                          .map(cls.new_allele)
            annotations['variant_allele'] = var_allele
            var_allele_freq = cls.dbsnp_allele_freq(annotations, variant_df)
            if var_allele_freq is not None:
                annotations['variant_allele_freq'] = var_allele_freq

        return annotations

    @staticmethod
    def cadd_genenames(cell):
        return [d['genename'] for d in cell
                if 'genename' in d]

    @staticmethod
    def listify(cell):
        if type(cell) != list:
            return [cell]

        return cell

    @staticmethod
    def dbsnp_allele_freq(annotations, variant_df):
        variant_allele_series = annotations['variant_allele']
        allele_freqs = variant_df['dbsnp.alleles']

        df = pd.concat([variant_allele_series, allele_freqs], axis=1)
        for hgvs in df.index.get_level_values('_id'):
            row = df.loc[pd.IndexSlice[:, hgvs], :].dropna(axis=1, how='all')
            allele = row['variant_allele'].iloc[0]

            if 'dbsnp.alleles' in row:
                dbsnp_alleles = row['dbsnp.alleles'].iloc[0]
                allele_info = pd.DataFrame(dbsnp_alleles).set_index('allele')

                if allele == 'del':  # This is returned by #new_allele() below
                    # The shorter allele is the one with the deletion.
                    # I can't decide if this is hackish or ok.
                    shorter_allele = allele_info.index.min()
                    freq = allele_info.loc[shorter_allele]['freq']
                else:
                    freq = allele_info.loc[allele]['freq']

                df.loc[pd.IndexSlice[:, hgvs], 'variant_allele_freq'] = freq

        if 'variant_allele_freq' not in df:
            return None

        return df['variant_allele_freq']

    @staticmethod
    def new_allele(hgvs):
        if 'del' in hgvs:
            return 'del'

        return hgvs.split('>')[-1]

    @staticmethod
    def snpeff_ann(cell):
        return list({(d['effect'], d['putative_impact'], d['gene_name'])
                     for d in cell})

    @staticmethod
    def info_present_from(label, df):
        return len(df.filter(regex=label).columns) > 0

    @staticmethod
    def _parse_grasp_publications_field(s):
        if 'grasp.publications' not in s:
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

