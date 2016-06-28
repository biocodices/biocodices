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
