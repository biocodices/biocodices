# Installation

Get these programs and install them:
* FastQC: http://www.bioinformatics.babraham.ac.uk/projects/download.html
* ea-utils: https://code.google.com/archive/p/ea-utils/downloads
* GatK: https://www.broadinstitute.org/gatk/download/
* Picard Tools: https://github.com/broadinstitute/picard/releases/tag/2.3.0
* Samtools: https://sourceforge.net/projects/samtools/files/
* Vcftools: http://vcftools.sourceforge.net/downloads.html
* Bedtools: https://github.com/arq5x/bedtools2/releases
* SmartPCA (EIGENSOFT): http://www.hsph.harvard.edu/alkes-price/software/
* Admixture: https://www.genetics.ucla.edu/software/admixture/download.html

```
# Ptyhon ternary for triangle plots
conda config --add channels conda-forge
conda install python-ternary

pip install myvariant
pip install multiqc
```

Create the settings file `~/.biocodices/executables.yml` with paths to every executable. This is mine, for instance:
```yaml
admixture: /home/juan/software/admixture_linux-1.3.0/admixture
smartpca: /home/juan/software/eigensoft6/src/eigensrc/smartpca
fastqc: /home/juan/software/FastQC/fastqc
fastq-mcf: /home/juan/software/ea-utils.1.1.2-537/fastq-mcf
bedtools: /home/juan/software/bedtools2/bin/bedtools
samtools: /home/juansoftware/samtools-1.3/samtools
gatk: java -jar /home/juan/software/GenomeAnalysisTK/GenomeAnalysisTK.jar
vcftools: /usr/local/bin/vcftools
picard-tools: java -jar /home/juan/software/picard-tools-2.2.4/picard.jar
```

Create a `~/.biocodices/resources.yml` file where you will specify a base
directory for some resources and the filenames they have inside that directory. Mine looks like this:
```yaml
base_dir: /home/juan/biocodices/resources
illumina_adapters_file: illumina_adapters.fasta
gwas_catalog_v1.0.1: &gwas_catalog_default gwas_catalog_v1.0-associations_e84_r2016-05-08.tsv
gwas_catalog: *gwas_catalog_default
clinvar:
    disease_names: clinvar_disease_names
    variants:
        GRCh37: clinvar_20160502_GRCh37.vcf
        GRCh38: clinvar_20160502_GRCh38.vcf
```

You need to download some resources from the web.

* Browse GATK bundle ftp servers to get the reference genome. For the GRCh37
    version, this was the command I run: `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta.gz`. Warning: the decompressed file will weight ~3Gb.
    - After getting the reference genome, unzip it (`gunzip -d human_g1k_v37.fasta.gz`) and run `bwa index -a bwtsw human_g1k_v37.fasta` to create the index files that bwa will need.

* A file with the adapters to trim from your reads.

# Recipes
## The Sequencing object and its Samples
```python
# Instatiate a Sequencing object with the root dir for its data.
# The directory should have a 'results' subdir with reads files.
sequencing = Sequencing('~/MyProject/NGS0001')
samples = sequencing.samples()  # => A list of sample objects
sample = samples[0]
sample.reads_filenames()  # => The forward and reverse reads files
```
## Samples have reads and can act on them
Each Sample object will take care of creating its own directory under results
(and under its sequencing root directory) to store the sample-specific files
there.
```python
for sample in sequencing.samples():
    sample.results_dir  # => Each sample has its own results directory
    sample.analyze_reads()  # => Generates HTML reports for each reads file
```
