import mygene


def annotate_gene(gene_symbol):
    """Query MyGene.info for a given gene symbol. E.g. 'PARK2'.
    Returns the json that the sites gives."""
    mg = mygene.MyGeneInfo()
    query_result = mg.query(gene_symbol)
    return query_result['hits']
