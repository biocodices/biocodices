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

    def protein_changes_by_rs(self):
        """Biocodices-parkinsonDB specific merging of tables."""
        df = self.table('mutations_report')
        df = df.loc[:, ['VariationName', 'GeneID','TypeVariation', 'GenomicChange', 'ProteinChange']]

        def parse_protein_change(change, ret):
            match = re.search(r'p.(\D{3})(\d+)(\D{3}|=)', change)
            if not match: return ''
            data = {'normal_aminoacid': match.group(1),
                    'protein_position': match.group(2),
                    'new_aminoacid': match.group(3)}
            return data[ret]

        df['Normal Aminoacid'] = df['ProteinChange'].map(lambda code: parse_protein_change(code, 'normal_aminoacid'))
        df['New Aminoacid'] = df['ProteinChange'].map(lambda code: parse_protein_change(code, 'new_aminoacid'))
        df['Protein Position'] = df['ProteinChange'].map(lambda code: parse_protein_change(code, 'protein_position'))

        protein_changes_by_rs = df.sort_values(by='GeneID').reset_index(drop=True)

        return protein_changes_by_rs
