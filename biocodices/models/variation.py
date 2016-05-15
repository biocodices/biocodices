from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float

Base = declarative_base()  # FIXME: where does this go??

class Variation(Base):
    __tablename__ = 'variations'
    id = Column(Integer, primary_key=True)
    name = Column(String(100), nullable=False)
    source = Column(String(50))
    hg19_start = Column(Integer, nullable=False)
    hg19_end = Column(Integer, nullable=False)
    chromosome = Column(String(50), nullable=False)
    alleles = Column(String(50), nullable=False)
    strand = Column(Integer)
    ancestral_allele = Column(Integer)
    maf_1kg = Column(Float)
    note = Column(String(2000))

    def __repr__(self):
        s = '<{}({}, chr)>'
        return s.format(self.__class__.__name__, self.name, self.source)
