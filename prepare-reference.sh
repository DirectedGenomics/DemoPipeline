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

usage() {
  cat << EOF
Usage: $0 <reference-fasta> <output-directory>

Prepares a reference sequence for use with the demonstration pipeline. It:
  - Normalizes the fasta file to ensure consistent line lengths
  - Generates a fasta index (.fai) file
  - Generates a sequence dictionary (.dict) file
  - Generates the bwa index files

The output directory should NOT be the directory containing the existing fasta.
After execution the output directory will contain a new fasta file and all
necessary files for the pipeline to run.

The first three steps may take a few minutes each.  Generation of the bwa
index can take multiple hours.
EOF
}

if [[ ! $# -eq 2 ]]; then usage; echo; echo "Error: must provide exactly two arguments."; exit 1; fi

FASTA=$1
OUT=$2

if [[ ! -e $FASTA ]]; then usage; echo; echo "Error: file $FASTA does not exist."; exit 1; fi

################################################################################
# Environment setup
################################################################################
P=`dirname $0`
source $P/common.sh
initialize

################################################################################
# Prepare the reference
################################################################################
FILENAME=$(basename $FASTA)
NEW_FASTA=$OUT/$FILENAME

execute "mkdir -p $OUT"

execute "java -Xmx4g -jar $picard NormalizeFasta I=$FASTA O=$NEW_FASTA LINE_LENGTH=80"

execute "samtools faidx $NEW_FASTA"

execute "java -Xmx4g -jar $picard CreateSequenceDictionary R=$NEW_FASTA"

execute "$bwa index $NEW_FASTA"

banner "Reference Preparation Complete"