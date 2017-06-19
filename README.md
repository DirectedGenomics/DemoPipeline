[![License](https://img.shields.io/badge/license-MIT-blue.svg)](./blob/master/LICENSE)
[![Language](https://img.shields.io/badge/language-bash-brightgreen.svg)](https://www.gnu.org/software/bash/)

# NEBNext Direct Demo Pipeline

This repository provides a minimal pipeline to process data generated with _NEBNext Direct_ kits.  It serves two purposes:

1. Provides a pipeline to generate BAMs and somatic variants in VCF format that can be used as-is in a non-production setting
2. Documents clearly a set of processing steps that can be transferred into any pipelining environment

## Pipeline Overview

The pipeline is implemented as a simple BASH script that uses the following open-source software:

* [bwa](https://github.com/lh3/bwa) - specifically `bwa mem` for alignment of reads to the genome
* [picard](https://broadinstitute.github.io/picard/) - used for various conversions, sorting, etc.
* [fgbio](https://github.com/fulcrumgenomics/fgbio) - used for generating consensus reads and filtering somatic variants
* [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) - to call somatic variants from the reads

The pipeline has the following general structure:

1. Generate an unmapped BAM from input fastq files, including UMIs
2. Map the raw reads to the genome and mark duplicates
3. Generate consensus reads from the raw reads
4. Re-map the consensus reads to the genome
5. Call and filter variants


The following files are created by the pipeline:

* `raw.unmapped.bam`: A BAM file of all reads prior to mapping
* `raw.mapped.bam`: A BAM file of reads post mapping with bwa
* `raw.deduped.bam`: A BAM file of all reads after duplicate marking
* `raw.deduped.alignment_summary_metrics.txt`: A text file containing summary metrics about the raw reads
* `raw.deduped.hs_metrics.txt`: A text file containing metrics about the target enrichment in the raw reads
* `grouped.bam`: A BAM file of reads grouped together by read positions and UMI
* `consensus.unmapped.bam`: A BAM file of consensus reads prior to mapping
* `consensus.mapped.bam`: A BAM file of consensus reads post mapping with bwa
* `consensus.filtered.bam`: A BAM file of consensus reads after filtering to reduce errors
* `consensus.filtered.alignment_summary_metrics.txt`: A text file containing summary metrics about the filtered consensus reads
* `consensus.filtered.hs_metrics.txt`: A text file containing metrics about the target enrichment in the filtered consensus reads
* `raw.vcf`: A VCF of variants called from the deduped raw reads
* `consensus.vcf`: A VCF of variants called from the filtered consensus reads

A production pipeline might choose to keep only the `raw.deduped.bam` and the `consensus.filtered.bam` BAM files and delete the remaining ones, but the intermediates are retained in the demonstration pipeline so that they may be examined.

## Running the Pipeline

### Prerequisites
In order to run the pipeline you will need:

1. A linux or mac computer with at least 12GB of RAM
2. [Java](https://java.com/en/download/manual.jsp) version 1.8 installed and available
3. [R](https://www.r-project.org/) version 3.0 or higher
4. A working PERL installation
5. A reference FASTA file with bwa index, fasta index and sequence dictionary (see `prepare-reference.sh` for help with this)
6. A BED file of target regions

### Preparing a Reference Sequence

Several files need to be generated from the reference FASTA file in order for the pipeline to function.  The `prepare-reference.sh` script can generate all the necessary files.  Alternatively you may manually prepare the reference using:

1. `samtools faidx ref.fa` to generate the FASTA index
2. `bwa index ref.fa` to generate the BWA index
3. `java -jar picard.jar CreateSequenceDictionary` to create the `.dict` or sequence dictionary file

### Downloading and Installing the Pipeline

The pipeline can be retrieved either by cloning this repository:

```
git clone https://github.com/DirectedGenomics/DemoPipeline.git
```

or by downloading one of the prepackaged releases from the [Releases](./releases) page and unzipping it.

### Executing the pipeline

The pipeline is executed by running the `pipeline.sh` script.  If run with no arguments, or incorrect arguments, it will print out information on usage and options.  An example invocation follows:

```bash
./pipeline.sh      \
  -r /refs/hg19.fa \
  -o /scratch/out  \
  -l targets.bed   \
  -s NA12878       \
  -t8              \
  r1.fq.gz r2.fq.gz i2.fq.gz
```

The pipeline requires _three_ fastq files as input and they must be in the expected order: read 1, followed by read 2, followed by the index read containing the UMIs (I2).

During execution the pipeline will emit log messages regarding it's progress and will also emit the full commands used for each step.  The results of the pipeline can be duplicated exactly by copying and pasting these commands, in order, into a script or terminal!
