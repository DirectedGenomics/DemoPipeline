#!/bin/bash

################################################################################
# MIT License
# 
# Copyright (c) 2017 Directed Genomics
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##################################################################################

set -euo pipefail

################################################################################
# Argument processing & setup
################################################################################

REF=""
OUT=""
BED=""
SAMPLE="testsample"
THREADS=4
MIN=1

usage() {
  cat << EOF
Usage: $0 [options] -r <reference-fasta> -l <target-bed> -o <directory> <r1.fq> <r2.fq> <i2.fq>

Analyzes a set of fastq files (optionally gzipped) that were generated using the
NEBNext Direct protocol.

OPTIONS:
   -r <ref.fasta>  The path to the reference fasta file
   -o <output dir> The name of the output directory to write to
   -l <target.bed> The path to the target regions BED file
   -s <sample-name> The name of the sequenced sample (default: $SAMPLE)
   -m <min-reads> The minimum number of reads to form a consensus read (default: $MIN)
   -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
EOF
}


while getopts “r:o:l:s:t:m:” OPTION; do
  case $OPTION in
    r) REF=$OPTARG;;
    o) OUT=$OPTARG;;
    l) BED=$OPTARG;;
    s) SAMPLE=$OPTARG;;
    m) MIN=$OPTARG;;
    t) THREADS=$OPTARG;;
    h) echo "Unknown option ${OPTION}"; usage; exit;;
    [?]) usage; exit;;
  esac
done

FQ1=${@:$OPTIND:1}
FQ2=${@:$OPTIND+1:1}
FQU=${@:$OPTIND+2:2}

if [[ -z "$REF" ]]; then usage; echo; echo "Error: reference fasta parameter is required."; exit 1; fi
if [[ -z "$OUT" ]]; then usage; echo; echo "Error: output directory parameter is required."; exit 1; fi
if [[ -z "$BED" ]]; then usage; echo; echo "Error: target bed file parameter is required."; exit 1; fi
if [[ -z "$FQ1" ]]; then usage; echo; echo "Error: three fastq files must be supplied."; exit 1; fi
if [[ -z "$FQ2" ]]; then usage; echo; echo "Error: three fastq files must be supplied."; exit 1; fi
if [[ -z "$FQU" ]]; then usage; echo; echo "Error: three fastq files must be supplied."; exit 1; fi

################################################################################
# Environment setup
################################################################################
P=`dirname $0`
source $P/common.sh
initialize

################################################################################
# Paths that are created/used in the pipeline
################################################################################
interval_list="$OUT/targets.interval_list"
raw_unmapped_bam="$OUT/raw.unmapped.bam"
raw_mapped_bam="$OUT/raw.mapped.bam"
raw_deduped_base="$OUT/raw.deduped"
raw_deduped_bam="${raw_deduped_base}.bam"
grouped_bam="$OUT/grouped.bam"
cons_unmapped_bam="$OUT/consensus.unmapped.bam"
cons_mapped_bam="$OUT/consensus.mapped.bam"
cons_filtered_base="$OUT/consensus.filtered"
cons_filtered_bam="${cons_filtered_base}.bam"

################################################################################
# Run the pipeline
################################################################################

execute "mkdir -p $OUT"

banner "Preparing unmapped BAM, mapping & deduping reads..."

execute "java -Xmx2g -jar $fgbio FastqToBam" \
            "--input $FQ1 $FQ2 $FQU" \
            "--read-structures +T +T +M" \
            "--output $raw_unmapped_bam" \
            "--sample $SAMPLE --library $SAMPLE" \
            "--sort"

execute "java -Xmx256m -jar $picard SamToFastq INPUT=$raw_unmapped_bam F=/dev/stdout QUIET=true INTERLEAVE=true" \
        " | $bwa mem -t $THREADS -p $REF /dev/stdin" \
        " | java -Xmx2g -jar $picard MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true QUIET=true" \
        "    UNMAPPED=$raw_unmapped_bam ALIGNED=/dev/stdin O=$raw_mapped_bam" \
        "    R=$REF " \
        "    ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN" \
        "    RV=cd RV=ce ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false"

execute "java -Xmx4g -jar $picard MarkDuplicates CREATE_INDEX=true" \
        "I=$raw_mapped_bam O=$raw_deduped_bam M=$OUT/duplicate_metrics.txt BARCODE_TAG=RX"

################################################################################
banner "Generating and re-aligning consensus reads..."

execute "java -Xmx4g -jar $fgbio GroupReadsByUmi" \
        "--input $raw_mapped_bam --output $grouped_bam" \
        "--family-size-histogram $OUT/tag_family_sizes.txt" \
        "--strategy=adjacency"

execute "java -Xmx4g -jar $fgbio CallMolecularConsensusReads" \
        "--input $grouped_bam --output $cons_unmapped_bam" \
        "--min-reads=$MIN --min-input-base-quality=20" \
        "--error-rate-pre-umi=50 --error-rate-post-umi=30"

execute "java -Xmx256m -jar $picard SamToFastq INPUT=$cons_unmapped_bam F=/dev/stdout QUIET=true INTERLEAVE=true" \
        " | $bwa mem -t $THREADS -p $REF /dev/stdin" \
        " | java -Xmx2g -jar $picard MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true QUIET=true" \
        "    UNMAPPED=$cons_unmapped_bam ALIGNED=/dev/stdin O=$cons_mapped_bam" \
        "    R=$REF " \
        "    ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN" \
        "    RV=cd RV=ce ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false"

if [ $MIN -eq 1 ];  then QUAL=30; else QUAL="44"; fi

execute "java -Xmx4g -jar $fgbio FilterConsensusReads" \
        "--input $cons_mapped_bam --output $cons_filtered_bam" \
        "--ref $REF" \
        "--min-reads=$MIN --min-base-quality=$QUAL"

################################################################################
banner "Generating Metrics"

execute "java -Xmx4g -jar $picard BedToIntervalList " \
        "I=$BED O=$interval_list SD=$REF"

execute "java -Xmx4g -jar $picard CollectAlignmentSummaryMetrics" \
        "I=$raw_deduped_bam O=${raw_deduped_base}.alignment_summary_metrics.txt R=$REF"

execute "java -Xmx4g -jar $picard CollectHsMetrics" \
        "I=$raw_deduped_bam O=${raw_deduped_base}.hs_metrics.txt" \
        "PER_TARGET_COVERAGE=${raw_deduped_base}.per_target_coverage.txt" \
        "R=$REF TI=$interval_list BI=$interval_list"

execute "java -Xmx4g -jar $picard CollectAlignmentSummaryMetrics" \
        "I=$cons_filtered_bam O=${cons_filtered_base}.alignment_summary_metrics.txt R=$REF"

execute "java -Xmx4g -jar $picard CollectHsMetrics" \
        "I=$cons_filtered_bam O=${cons_filtered_base}.hs_metrics.txt" \
        "PER_TARGET_COVERAGE=${cons_filtered_base}.per_target_coverage.txt" \
        "R=$REF TI=$interval_list BI=$interval_list"

################################################################################
for datatype in "raw" "consensus"; do
    banner "Calling and filtering variants in $datatype reads..."

    if [ $datatype == "raw" ]; then 
        bam=$raw_deduped_bam
        vcf=$OUT/raw.vcf
    else
        bam=$cons_filtered_bam; 
        vcf=$OUT/consensus.vcf
    fi

    execute "$vddir/bin/VarDict" \
            "-G $REF" \
            "-N $SAMPLE -b $bam" \
            "-z 1 -c 1 -S 2 -E 3 -g 4 -F 0x700 -f 0.005" \
            "-r 3 -q 30 -th $THREADS $BED" \
            "| awk '{if (\$6 != \$7) print}'" \
            "| $vddir/bin/teststrandbias.R" \
            "| $vddir/bin/var2vcf_valid.pl -N $SAMPLE -E -f 0.01" \
            "> $OUT/tmp.vcf"
    
    fasta_ext=${REF##*.}
    fasta_dir=$(dirname $REF)
    dict="${fasta_dir}/$(basename $REF $fasta_ext)dict"
    
    execute "java -Xmx2g -jar $picard SortVcf I=$OUT/tmp.vcf O=$OUT/sorted.vcf CREATE_INDEX=true SD=$dict"
    
    execute "java -Xmx2g -jar $fgbio FilterSomaticVcf" \
            "--input $OUT/sorted.vcf" \
            "--output $vcf" \
            "--bam $bam" \
            "--min-mapping-quality 1" \
            "--end-repair-distance 15" \
            "--end-repair-p-value 0.05"

    execute "rm -f $OUT/tmp.vcf $OUT/sorted.vcf"
done

banner "Completed."
