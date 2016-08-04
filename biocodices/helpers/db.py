import re
import pandas as pd
from sqlalchemy import create_engine

from biocodices.helpers import Config


class DB:
    def __init__(self, database):
        """
        Create a connection to a db by name. Use the DB instance to query it.
        """
        creds = Config('mysql_credentials')
        creds.update({'database': database})
        url = 'mysql+pymysql://{user}:{pass}@{host}/{database}'.format(**creds)
        engine = create_engine(url)

        self.conn = engine.connect()

        # Store the tables in memory and close the DB connection
        for table_name in self.tables:
            setattr(self, table_name, self.table(table_name))

        if database == 'parkinsonDB':
            self.enp1_variants = self._merge_enp1_tables()

    def query(self, query_string, as_list=False):
        result = self.conn.execute(query_string)
        if as_list:
            result = list(result)
        return result

    def table(self, table_name):
        """Retrieve by name a database table as a pandas DataFrame."""
        return pd.read_sql_table(table_name, self.conn)

    @property
    def tables(self):
        """Show table names from the current database."""
        return [t[0] for t in self.query('show tables;', as_list=True)]

    def _merge_enp1_tables(self):
        """Merge the ENPv1 parkinsonDB tables to get the most relevant data
        for each panel marker. Returns a big DataFrame."""
        df = self.protein_changes_by_rs()
        original_columns = list(df.columns)

        mutrep_columns = ['VariationName', 'RISK_ALLELE']
        df = pd.merge(df, self.mutations_report[mutrep_columns],
                      on='VariationName').drop_duplicates()

        muteff_columns = ['VariationName', 'Polyphen_prediction',
                          'Sift_prediction']
        df = pd.merge(df, self.mutations_effects[muteff_columns],
                      on='VariationName').drop_duplicates()

        cit_columns = ['VariationName', 'Review', 'Note']
        df = pd.merge(df, self.variationscitations[cit_columns],
                      on='VariationName').drop_duplicates()

        # Set column order following the order of DFs merged
        columns = original_columns
        columns += [col for col in mutrep_columns if col != 'VariationName']
        columns += [col for col in muteff_columns if col != 'VariationName']
        columns += [col for col in cit_columns if col != 'VariationName']
        df = df[columns]

        new_df = pd.DataFrame({})

        for variant, variant_df in df.groupby('VariationName'):
            new_row = variant_df.reset_index(drop=True)
            new_row = new_row.loc[0, 'VariationName':'RISK_ALLELE'].copy()

            polyphen_predictions = variant_df['Polyphen_prediction'].unique()
            polyphen_predictions = [pred for pred in polyphen_predictions
                                    if pred and pred != '']
            new_row['Polyphen_prediction'] = ','.join(polyphen_predictions)

            sift_predictions = variant_df['Sift_prediction'].unique()
            sift_predictions = [pred for pred in sift_predictions
                                if pred and pred != '']
            new_row['Sift_prediction'] = ','.join(sift_predictions)

            reviews = variant_df[['Review', 'Note']].drop_duplicates()['Review']
            reviews = [review.strip() for review in reviews if review != '']
            new_row['Review'] = ','.join(reviews)

            new_df = new_df.append(new_row, ignore_index=True)

        column_order = [col for col in df.columns if col != 'Note']
        new_df = new_df.sort_values(by='GeneID')

        return new_df[column_order].set_index('VariationName')

    def protein_changes_by_rs(self):
        """Biocodices-parkinsonDB specific merging of tables."""
        df = self.mutations_report
        df = df.loc[:, ['VariationName', 'GeneID','TypeVariation', 'GenomicChange', 'ProteinChange']]

        def parse_protein_change(change, ret):
            match = re.search(r'p.(\D{3})(\d+)(\D{3}|=)', change)
            if not match: return ''

            data = {'normal_aminoacid': match.group(1),
                    'protein_position': match.group(2),
                    'new_aminoacid': match.group(3).replace('Ter', 'STOP')}
            return data[ret]

        df['Normal Aminoacid'] = df['ProteinChange'].map(lambda code: parse_protein_change(code, 'normal_aminoacid'))
        df['New Aminoacid'] = df['ProteinChange'].map(lambda code: parse_protein_change(code, 'new_aminoacid'))
        df['Protein Position'] = df['ProteinChange'].map(lambda code: parse_protein_change(code, 'protein_position'))

        protein_changes_by_rs = df.sort_values(by='GeneID').reset_index(drop=True)
        ret = protein_changes_by_rs.set_index('VariationName').join(self.alleles_by_rs())

        return ret.reset_index()

    def alleles_by_rs(self):
        """Biocodices-parkinsonDB specific merging of tables."""
        alleles = self.variations.set_index('VariationName')
        alleles['Alleles'] = alleles['Alleles'].str.split('/')
        alleles['ancestral_allele'] = alleles['Alleles'].map(lambda l: l[0])
        alleles['new_alleles'] = alleles['Alleles'].map(lambda l: ','.join(l[1:]))

        return alleles[['ancestral_allele', 'new_alleles']]
