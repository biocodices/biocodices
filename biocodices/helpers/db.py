import pandas as pd
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from sqlalchemy import create_engine

from biocodices.helpers import Config


class DB:
    def __init__(self, database):
        """
        Create a connection to a db by name. Use the DB instance to query it.
        """
        Base = automap_base()
        creds = Config('mysql_credentials')
        creds.update({'database': database})
        url = 'mysql+pymysql://{user}:{pass}@{host}/{database}'.format(**creds)
        engine = create_engine(url)

        # reflect the existing tables
        Base.prepare(engine, reflect=True)

        # Hardcoded from our DB. This should go in a settings file.
        if database == 'parkinsonDB':
            self.variation = Base.classes.variations
            self.variation_citation = Base.classes.variationscitations
            self.gene = Base.classes.gene
        if database == 'biocodices':
            self.sample = Base.classes.MUESTRA
            self.ngs_sample = Base.classes.NGS_MUESTRAS

        self.session = Session(engine)
        self.conn = engine.connect()

    def query(self, query_string, as_list=False):
        result = self.conn.execute(query_string)
        if as_list:
            result = list(result)
        return result

    def table(self, table_name):
        """Retrieve a database table as a pandas DataFrame."""
        return pd.read_sql_table(table_name, self.conn)

    def table_names(self):
        return [t[0] for t in self.query('show tables;', as_list=True)]
