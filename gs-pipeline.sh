#!/bin/bash

################################################################################
# MIT License
#
# Copyright (c) 2020 Directed Genomics
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
GATK_VERSION="4.1.6.0"
GATK_URL="https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip"

usage() {
  cat << EOF
Usage: $0 [options] -r <reference-fasta> -l <target-bed> -o <directory> <r.fq> <u.fq>

Analyzes a set of fastq files (optionally gzipped) that were generated using the
NEBNext Direct GS protocol.

OPTIONS:
   -r <ref.fasta>  The path to the reference fasta file
   -o <output dir> The name of the output directory to write to
   -l <target.bed> The path to the target regions BED file
   -s <sample-name> The name of the sequenced sample (default: $SAMPLE)
   -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
EOF
}


while getopts “r:o:l:s:t:m:” OPTION; do
  case $OPTION in
    r) REF=$OPTARG;;
    o) OUT=$OPTARG;;
    l) BED=$OPTARG;;
    s) SAMPLE=$OPTARG;;
    t) THREADS=$OPTARG;;
    h) echo "Unknown option ${OPTION}"; usage; exit;;
    [?]) usage; exit;;
  esac
done

FQ1=${@:$OPTIND:1}
FQU=${@:$OPTIND+1:1}

if [[ -z "$REF" ]]; then usage; echo; echo "Error: reference fasta parameter is required."; exit 1; fi
if [[ -z "$OUT" ]]; then usage; echo; echo "Error: output directory parameter is required."; exit 1; fi
if [[ -z "$BED" ]]; then usage; echo; echo "Error: target bed file parameter is required."; exit 1; fi
if [[ -z "$FQ1" ]]; then usage; echo; echo "Error: read fastq file must be supplied."; exit 1; fi
if [[ -z "$FQU" ]]; then usage; echo; echo "Error: UMI fastq file must be supplied."; exit 1; fi

################################################################################
# Environment setup
################################################################################
P=$(dirname $0)
source $P/common.sh
initialize

################################################################################
# Paths that are created/used in the pipeline
################################################################################
umi_out="$OUT/UMI"
interval_list="$OUT/targets.interval_list"
raw_unmapped_bam="$OUT/$SAMPLE.unmapped.bam"
raw_mapped_bam="$OUT/$SAMPLE.mapped.bam"
raw_deduped_base="$umi_out/$SAMPLE.deduped"
raw_deduped_bam="${raw_deduped_base}.bam"
grouped_bam="$OUT/grouped.bam"


################################################################################
# Run the pipeline
################################################################################

if [ ! -f $P/bin/gatk.jar ]; then
    banner "Downloading GATK4..."
    if $(which wget > /dev/null); then
         wget -qO $P/gatk.zip $GATK_URL
    elif $(which curl > /dev/null); then
        curl -Lso $P/gatk.zip $GATK_URL
    else
        fail "wget or curl must be installed and available in order to download GATK."
    fi

    unzip -p $P/gatk.zip gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-local.jar > $P/bin/gatk.jar
    rm $P/gatk.zip
fi

execute "mkdir -p $OUT"
execute "mkdir -p $umi_out"

banner "Preparing unmapped BAM, mapping & deduping reads..."

execute "java -Xmx2g -jar $fgbio FastqToBam" \
            "--input $FQ1 $FQU" \
            "--read-structures +T +M" \
            "--output $raw_unmapped_bam" \
            "--sample $SAMPLE --library $SAMPLE" \
            "--sort"

execute "java -Xmx256m -jar $picard SamToFastq INPUT=$raw_unmapped_bam F=/dev/stdout QUIET=true INTERLEAVE=true" \
        " | $bwa mem -t $THREADS -p $REF /dev/stdin" \
        " | java -Xmx2g -jar $picard MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true QUIET=true" \
        "    UNMAPPED=$raw_unmapped_bam ALIGNED=/dev/stdin O=/dev/stdout" \
        "    R=$REF " \
        "    ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN" \
        "    RV=cd RV=ce ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false" \
        " |  java -Xmx2g -jar $fgbio FilterBam" \
        " --input /dev/stdin --output $raw_mapped_bam --remove-duplicates false --min-map-q 0" \
        " --remove-secondary-alignments false"


for datatype in "raw" "dedup"; do
    banner "Calling and filtering variants in $datatype reads..."
    in_bam=$raw_mapped_bam

    if [ $datatype == "raw" ]; then
        bam=$raw_mapped_bam
        out_dir=$OUT
        germline_vcf=$OUT/$SAMPLE.vcf.gz
    else
        germline_vcf=$umi_out/$SAMPLE.vcf.gz
        bam=$raw_deduped_bam
        out_dir=$umi_out

        execute "java -Xmx4g -jar $picard MarkDuplicates CREATE_INDEX=true" \
                "I=$in_bam O=$bam M=$umi_out/duplicate_metrics.txt BARCODE_TAG=RX"
    fi

    out_prefix=$out_dir/$SAMPLE

    banner "Generating Metrics"

    execute "java -Xmx4g -jar $picard BedToIntervalList " \
            "I=$BED O=$interval_list SD=$REF"

    execute "java -Xmx4g -jar $picard CollectAlignmentSummaryMetrics" \
            "I=$bam O=${out_prefix}.alignment_summary_metrics.txt R=$REF"

    execute "java -Xmx4g -jar $picard CollectHsMetrics" \
            "I=$bam O=${out_prefix}.hs_metrics.txt" \
             "PER_TARGET_COVERAGE=${out_prefix}.per_target_coverage.txt" \
            "R=$REF TI=$interval_list BI=$interval_list"

    # Germline Calling
    execute "java -Xmx4g -jar $gatk HaplotypeCaller" \
            "--reference $REF" \
            "--intervals $BED" \
            "--minPruning 3" \
            "--maxNumHaplotypesInPopulation 200" \
            "--emitRefConfidence GVCF" \
            "--max_alternate_alleles 3" \
            "--contamination_fraction_to_filter 0.0" \
            "--input $bam" \
            "--output $out_dir/tmp.g.vcf.gz"

    execute "$tabix $out_dir/tmp.g.vcf.gz"

    execute "java -Xmx4096m -jar $gatk GenotypeGVCFs" \
            "--reference $REF" \
            "--intervals $BED" \
            "--variant $out_dir/tmp.g.vcf.gz" \
            "--output $out_dir/tmp.germline.unfiltered.vcf.gz"

    execute "java -Xmx4096m -jar $picard FilterVcf" \
             "VALIDATION_STRINGENCY=SILENT" \
             "CREATE_INDEX=true" \
             "I=$out_dir/tmp.germline.unfiltered.vcf.gz" \
             "O=$germline_vcf" \
             "MIN_AB=0.2 MIN_DP=0 MIN_GQ=20 MAX_FS=50.0 MIN_QD=6.0"

done

banner "Completed."
