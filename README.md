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

# Client for MyVariant.info
pip install myvariant
```

Create the settings file `~/.biocodices/executables.yml` with paths to every executable. This is mine, for instance:
```
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

Create a `~/.biocodices/resources.yml` file. Mine looks like this:
```

```




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
