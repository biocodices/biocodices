class Variation(Base):
    __tablename__ = 'variations'
    id = Column(Integer, primary_key=True)
    name = Column(String(100), nullable=False)
    source = Column(String(50))
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    chromosome = Column(String(50), nullable=False)
    alleles = Column(String(50), nullable=False)
    strand = Column(Integer)
    ancestral_allele = Column(Integer)
    maf_1kg = Column(Float)
    note = Column(String(2000))

