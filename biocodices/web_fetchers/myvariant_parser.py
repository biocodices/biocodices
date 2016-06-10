import pandas as pd


class MyvariantParser:
    @classmethod
    def publications(cls, myvariant_series):
        publications = []
        publications += cls._parse_grasp_publications_field(myvariant_series)
        return publications

    @classmethod
    def parse_annotations(cls, variant_df):
        if variant_df.empty:
            return pd.DataFrame({})

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
            # NOTE: cadd.gene.genename and cadd.gene are never together, and
            # whenever info from one of those fields is present, the other is
            # absent. That's why I merge them here in the same key
            # 'cadd.gene.genename', to deal with one column for this CADD datum.
            gathered_info['cadd.gene.genename'] = \
                variant_df['cadd.gene'].dropna().map(cls.cadd_genenames)

        if 'dbsnp.chrom' in variant_df:
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

        if '_id' in variant_df.index.names:
            var_allele = annotations.index.get_level_values('_id') \
                                          .map(cls.new_allele)
            annotations['variant_allele'] = var_allele

        if 'dbsnp.alleles' in variant_df:
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

                if allele_info.empty:
                    # This happens when there's info about different alleles
                    # But NOT a frequency associated to any of them in the
                    # dictionary 'dbsnp_alleles'.
                    return None

                if allele == 'del':  # This is returned by #new_allele() below
                    # The shorter allele is the one with the deletion.
                    # I can't decide if this is hackish or ok.
                    shorter_allele = allele_info.index.min()
                    freq = allele_info.loc[shorter_allele]['freq']
                elif allele == 'ins':
                    longer_allele = allele_info.index.max()
                    freq = allele_info.loc[longer_allele]['freq']
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

        if 'ins' in hgvs:
            return 'ins'

        return hgvs.split('>')[-1]

    @staticmethod
    def snpeff_ann(cell):
        return list({(d['effect'], d['putative_impact'], d['gene_name'])
                     for d in cell})

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

