from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from biocodices.helpers.config import Config

## Don't commit this
def create_dbsession(credentials=None):
    credentials = credentials or Config('mysql_credentials')
    conn_str = 'mysql+pymysql://{user}:{pass}@{host}/parkinsonDB'
    engine = create_engine(conn_str.format(**credentials))
    Session = sessionmaker(bind=engine)
    session = Session()

    return session
