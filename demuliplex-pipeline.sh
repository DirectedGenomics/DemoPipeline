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

RUN=""
OUT=""
SAMPLEFILE=""
READSTRUCT=""
THREADS=4
MERGE=0
MIN=1

usage() {
  cat << EOF
Usage: $0 [options] -r <run dir> -s <sample file> -o <directory>

Demultiplexes a sequencing run of libraries created from the NEBNext Direct GS protocol.  Sample and barcode
information to be used for demultiplexing is provided by a text file listing the samples and barcodes. See
the sample-file-template for the required format. All columns must be filled in for each sample present, regardless
if the sample has zero, one, or two barcodes. For cases where no barcode is present put an "N" in the appropriate
barcode column.

OPTIONS:
   -r <run dir>  The path to the sequencing run directory
   -o <output dir> The name of the output directory to write to
   -s <sample file> The file containing names of the samples and barcodes to demultiplex
   -p <read structure> The read structure used to translate the bases in a sequencing run into reads
   -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
   -m <merge lanes> Merge lanes (default: $MERGE)
EOF
}


while getopts “r:o:s:p:t:m:” OPTION; do
  case $OPTION in
    r) RUN=$OPTARG;;
    o) OUT=$OPTARG;;
    s) SAMPLEFILE=$OPTARG;;
    p) READSTRUCT=$OPTARG;;
    t) THREADS=$OPTARG;;
    m) MERGE=$OPTARG;;
    h) echo "Unknown option ${OPTION}"; usage; exit;;
    [?]) usage; exit;;
  esac
done


if [[ -z "$RUN" ]]; then usage; echo; echo "Error: run directory parameter is required."; exit 1; fi
if [[ -z "$OUT" ]]; then usage; echo; echo "Error: output directory parameter is required."; exit 1; fi
if [[ -z "$SAMPLEFILE" ]]; then usage; echo; echo "Error: sample file parameter is required."; exit 1; fi
if [[ -z "$READSTRUCT" ]]; then usage; echo; echo "Error: read structure parameter is required."; exit 1; fi


################################################################################
# Environment setup
################################################################################
P=$(dirname $0)
source $P/common.sh
initialize

execute "mkdir -p $OUT"

num_samples=$(egrep [AGCTN] $SAMPLEFILE | wc -l |  awk '{print $1}')
num_lanes=$(xpath $RUN/RunInfo.xml "string(RunInfo/Run/FlowcellLayout/@LaneCount)" 2> /dev/null)
flowcell=$(xpath $RUN/RunInfo.xml "RunInfo/Run/Flowcell/text()" 2> /dev/null)
num_barcodes=$(echo $READSTRUCT | tr -cd "B" | wc -c)

#write lib params file header
header="OUTPUT_PREFIX\tSAMPLE_ALIAS\tLIBRARY_NAME"
extr_header="OUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME"
if [ $num_barcodes -ge 1 ]; then
    header="$header\tBARCODE_1"
    extr_header="$extr_header\tbarcode_sequence_1"
fi
if [ $num_barcodes -ge 2 ]; then
    header="$header\tBARCODE_2"
    extr_header="$extr_header\tbarcode_sequence_2"
fi

for lane in $(seq 1 $num_lanes); do
    banner "Demultiplexing lane $lane ..."

    lib_params_file=$OUT/library_params.l"$lane".txt
    lib_extract_file=$OUT/library_params_extract.l"$lane".txt

    echo -e "$header" > $lib_params_file
    echo -e "$extr_header" > $lib_extract_file

    {
    read
    while read -r sample barcode1 barcode2;
    do
        values="${OUT}/${sample}.l${lane}\t${sample}\t${sample}"
        if [ $num_barcodes -ge 1 ]; then
          values="$values\t$barcode1"
        fi
        if [ $num_barcodes -ge 2 ]; then
          values="$values\t$barcode2"
        fi
        echo -e $values >> $lib_params_file
        echo -e $values >> $lib_extract_file
    done
    } < $SAMPLEFILE

    if [ $num_barcodes -gt 0 ]; then
        execute "java -Xmx4g -jar $picard ExtractIlluminaBarcodes" \
                "VALIDATION_STRINGENCY=SILENT BASECALLS_DIR=$RUN/Data/Intensities/BaseCalls/" \
                "OUTPUT_DIR=$RUN/Data/Intensities/BaseCalls/ READ_STRUCTURE=$READSTRUCT BARCODE_FILE=$lib_extract_file" \
                "METRICS_FILE=$RUN/Data/Intensities/BaseCalls/barcode_counts.lane-$lane.metrics.txt" \
                "COMPRESS_OUTPUTS=true NUM_PROCESSORS=$THREADS MAX_MISMATCHES=2 LANE=$lane"
    fi

    execute "java -Xmx12g -jar $picard IlluminaBasecallsToFastq" \
            "VALIDATION_STRINGENCY=SILENT BASECALLS_DIR=$RUN/Data/Intensities/BaseCalls/" \
            "READ_STRUCTURE=$READSTRUCT MULTIPLEX_PARAMS=$lib_params_file" \
            "INCLUDE_NON_PF_READS=false IGNORE_UNEXPECTED_BARCODES=true LANE=$lane RUN_BARCODE=$flowcell" \
            "NUM_PROCESSORS=$THREADS READ_NAME_FORMAT=ILLUMINA"
done

banner "Completed."