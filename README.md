# Installation

## Python and some libraries

(TODO: This section needs more detail, conda, etc.)

```
# Ptyhon ternary for triangle plots
conda config --add channels conda-forge
conda install python-ternary

pip install myvariant
pip install multiqc
```

## Software and resources

Get these programs and install them:

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/download.html)
* [ea-utils](https://code.google.com/archive/p/ea-utils/downloads)
* [GatK](https://www.broadinstitute.org/gatk/download/)
* [Picard Tools](https://github.com/broadinstitute/picard/releases/tag/2.3.0)
* [Samtools](https://sourceforge.net/projects/samtools/files/)
* [Vcftools](http://vcftools.sourceforge.net/downloads.html)
* [Bedtools](https://github.com/arq5x/bedtools2/releases)
* [SmartPCA (EIGENSOFT)](http://www.hsph.harvard.edu/alkes-price/software/)
* [Admixture](https://www.genetics.ucla.edu/software/admixture/download.html)

You also need to download some resources from the web:

Browse GATK bundle ftp servers to get the reference genome. For the GRCh37 version, I ran this command: `wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta.gz`. **Warning**: the decompressed file will weight ~3Gb.

Also from GATK servers, get files for known indels:
    - 1000 Genomes indels
    - Mills and 1000 Genomes Gold Standard
    - All known variants in GRCh37

```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase1.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf.gz
```

Unzip all the GATK bundle files with a `gunzip <filename>` command.

Prepare the reference genome fasta for BWA and GATK use, as detailed [here](http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference). You need the software installed in the first step for this:
```bash
bwa index -a bwtsw human_g1k_v37.fasta
samtools faidx human_g1k_v37.fasta
java -jar picard.jar CreateSequenceDictionary \
    REFERENCE=human_g1k_v37.fasta \
    OUTPUT=human_g1k_v37.dict
```

Finally, create a fasta file with the adapters to trim from your reads. I got the sequences from [Illumina's TruSeq documentation](http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf). Make sure you create a valid fasta file ('>ID' in one line, sequence in the next one, repeat). The IDs don't really matter, the sequences will be read.

# Settings before running `biocodices`

`biocodices` will search for its config files in a `~/.biocodices` directory,
so create it with `mkdir -p ~/.biocodices`. You will also need an organized
resource folder, so I recommend creating a `~/biocodices/resources` directory:
`mkdir -p ~/biocodices/resources`. Put the resources you downloaded from the
web in there.

Create the settings file `~/.biocodices/executables.yml` with paths to every executable from the software you downloaded earlier. Avoid tildes (`~`), write whole paths. This is my file, for instance:

```yaml
admixture: /home/juan/software/admixture_linux-1.3.0/admixture
smartpca: /home/juan/software/eigensoft6/src/eigensrc/smartpca
fastqc: /home/juan/software/FastQC/fastqc
fastq-mcf: /home/juan/software/ea-utils.1.1.2-537/fastq-mcf
bedtools: /home/juan/software/bedtools2/bin/bedtools
samtools: /home/juansoftware/samtools-1.3/samtools
GATK: java -jar /home/juan/software/GenomeAnalysisTK/GenomeAnalysisTK.jar
vcftools: /usr/local/bin/vcftools
picard-tools: java -jar /home/juan/software/picard-tools-2.2.4/picard.jar
```

Create a `~/.biocodices/resources.yml` file where you will specify a base
directory for the resources and the filenames they have inside that directory.
Mine looks like this:

```yaml
base_dir: /home/juan/biocodices/resources
illumina_adapters_file: illumina_adapters.fasta
panel_amplicons:
    ENPv1: ENPv1_amplicons_sorted.bed
reference_genome_hg19: &reference_genome_default human_g1k_v37.fasta
reference_genome: *reference_genome_default
1000G_indels: 1000G_phase1.indels.b37.vcf
mills_indels: Mills_and_1000G_gold_standard.indels.b37.vcf
gwas_catalog_v1.0.1: &gwas_catalog_default gwas_catalog_v1.0-associations_e84_r2016-05-08.tsv
gwas_catalog: *gwas_catalog_default
clinvar:
    disease_names: clinvar_disease_names
    variants:
        GRCh37: clinvar_20160502_GRCh37.vcf
        GRCh38: clinvar_20160502_GRCh38.vcf
```

# Recipes

## Feed `biocodices` with the input `fastq` files

Create a directory for a given sequencer run and create a `data` subdirectory
in it. Put the forward and reverse `fastq` files from all samples in `data`.
`biocodices` expects them to be named following this pattern:
`<sample_ID>.R1.fastq` for the forward reads, and `<sample_ID>.R2.fastq` for
the reverse reads of the same sample.

After putting the input `fastq` files in there, you're good to go!

## The Sequencing object and its Samples

```python
# Instatiate a Sequencing object with the root dir of a sequencer run.
# The directory should have a 'data' subdir with the fastq files.
sequencing = Sequencing('~/MyProject/NGS0001')
samples = sequencing.samples()  # => A list of sample objects
sample = samples[0]
sample.reads_filenames()  # => The forward and reverse reads files
```

## Variant calling

Each Sample object will take care of creating its own results directory
(named after the sample ID under the sequencing results directory).

You can call variants for each sample easily:

```python
for sample in sequencing.samples():
    sample.call_variants()
```

The process will take some times (i.e. 5 mins per sample, YMMV) and each step
of the process will be logged to a different log file under the sample's result
directory. To check where that is, you can ask: `sample.results_dir`, but in
general results and logs will be stored in `<sequencing_dir>/results/<sample_dir>/`.
