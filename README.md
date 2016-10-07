# Biocodices

## Annotation utilities

```python
from biocodices.annotation import Ensembl
from biocodices.annotation import DbSNP
from biocodices.annotation import MyVariantAnnotator

rs_list = ['rs21229385', 'rs21229338']

ensembl = Ensembl()
annotations = ensembl.annotate(rs_list)

dbsnp = DbSNP()
annotations = dbsnp.annotate(rs_list)

myvariant = MyVariantAnnotator()
annotations = myvariant.annotate(rs_list)
```

Clinvar annotator needs Clinvar IDs to query, so a useful move might involve
annotating a list of rs IDs with dbsnp, Ensembl or MyVariant and checking
if they have an associated Clinvar ID.

```python
from biocodices.annotation import Clinvar

clinvar = Clinvar()
annotations = clinvar.annotate(clinvar_ids)
```
