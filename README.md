# Variant calling

## Raw reads (fastq):

* Fastq: Quitar las secuencias repetidas, adaptadores...
  Recalcular la calidad por base.
* Fastq: Alinear con el genoma de referencia: BAM.
* Filtrar por calidad? BAM.
* Realign indels: BAM.
* Recalibrate bases after realignment: BAM.

## Analysis-ready reads (BAM):
* Variant calling: vcf
* Joint genotyping
* Variant Filtering

# PANGEA
## Analysis-ready variants (vcf):
* Variant annotation:
    - SnpEff
    - VEP
    - COSMIC
    - ...

* Report

# Protocols
## Design the case-control study

1. Specify case definition
2. Determine if the disease is heritable: if it's < 20% and the disease is common
   it's likely that very large samples will be required.
3. Is a population-based approach the right choice? It is if the disease of interest
   could reasonably include one or more *common underlying polymorphisms*.
4. Control selection: 'classical epidemiological designs', and collect relevant
   phenotypic and covariate info for the controls.
5. Determine the required sample size:
    i. Determine a minimum odds ratio of the disease allele to be detected in the study.
       When replicating a previous association, use a smaller odds ratio than that
       from the hypothesis-generating study.

## Marker selection

1. Visualize genomic information for a candidate gene to select tagSNPs.
   UCSC Browser. You can search by position or identifier(s).
    i.  Choose the RefSeq gene that covers the largest area of the gene of interest.
        Keep track of the version and date of access of the browser! Data will change.
    ii. Visualize the LD structure around the candidate gene. Expand the region
        of interest either because LD continues to the sides or just by 10kb to
        include nearby regulatory regions.
    iii. Retrieve HapMap SNPs in the region.
    iii. Create a list of functional SNPs to be forced into the selection.
         Use UCSC functional annotation.


       # - Check the startpos and endpos of the gene with UCSC table browser
       # - Use Haploview to download a +- 50kb region of LD and check LD there?
       # - With the new coordinates, download SNP genotype data from UCSC?
       # - Use UCSC table browser Variation group to dload SNPs in a region

## Statistical analysis of case-control studies
1. From the .map/.ped files, run the tests with PLINK
    * plink --file <label> --assoc --out <out_label>
      The genotypes must already be filtered by individuals with missing data rate of < 20%,
      SNPs with missing rate of > 5%, MAF <1%, HWE P-value < 1e-4
    * plink --file <label> --model --out <out_label>
      Runs all association tests: TREND, ALLELIC, GENO, RECESSIVE, DOMINANT
2. Visualize the result of the tests: manhattan plot, LD plot with Haploview
3. Multiple testing adjustment:
    * Bonferroni (alpha / n)
    * Holm, Sidak and FDR: plink --file <label> --adjust --model-{trend,rec,dom} --out <outlabel>
    * Permutation testing: --model --mperm 1000 --model-{trend,rec,dom} --out ...
4. Check for population stratification
