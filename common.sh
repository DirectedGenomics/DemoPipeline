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

# Initializes the pipeline environment so that scripts can run as expected
function initialize() {
    IFS=$'\n\t'
    detect_platform
    bwa=$P/bin/$PLATFORM/bwa
    tabix=$P/bin/$PLATFORM/tabix
    picard=$P/bin/picard.jar
    fgbio=$P/bin/fgbio.jar
    gatk=$P/bin/gatk-package-4.beta.1-local.jar
    vddir=$P/bin/VarDictJava
    check_dependencies
}

# Function to check and ensure that the dependencies necessary to run the pipeline
# are available
function check_dependencies() {
    # Check that an appropriate version of Java is available
    if ! type -p java > /dev/null ; then
        fail "Could not find java executble. Please ensure Java 1.8 or above is on the path." \
             "To download the latest version go to: https://java.com/en/download/manual.jsp"
    fi
    
    java_version=$(java -version 2>&1 | fgrep version | cut -d\  -f3 | tr -d '"' | cut -d. -f1-2)
    major=$(echo $java_version | cut -d. -f1)
    minor=$(echo $java_version | cut -d. -f2)
    
    if [ $major -lt 2 ] && [ $minor -lt 8 ]; then
        fail "Detected java version $java_version on the path. Java version 1.8+ is required" \
             "to run this pipeline. To download the latest version go to: " \
             "https://java.com/en/download/manual.jsp"
    fi
    
    # TODO: CHECK THAT R is available
}

# Function to detect that platform being executed on and set PLATFORM
function detect_platform() {
  if   $(uname -a | fgrep -i darwin > /dev/null); then PLATFORM="mac"
  elif $(uname -a | fgrep -i linux  > /dev/null); then PLATFORM="linux"
  else fail "Could not detect supported operating system."
  fi
}

#
# Function used to exit after printing a large error message
function fail() {
    banner $*
    exit 1
}

function execute() {
    log "--------------------------------------------------------------------------------"
    log "- Executing: "
    for i in `seq 1 $#`; do
        if [ $i -eq $# ]; then lineterm=""; else lineterm=" \\"; fi
        if [ $i -eq 1 ];  then prefix=""; else prefix="    "; fi
        log "- ${prefix}${!i}${lineterm}"
    done    
    log "--------------------------------------------------------------------------------"
    
    OLD=$IFS
    IFS=" "
    command="$@"
    eval $command
    IFS=$OLD
}

# Logs a message to the console along with the date an time
function log() {
    echo [`date`] $*
}

# Simple function to take arguments and return them as a file-system-safe string
function make_fss() {
    echo $* | tr '!$#()[]' '-'
}

# A short function for echoing a string to the screen in a banner
function banner() {
    echo
    echo "################################################################################"
    for line in "$@"; do
        echo "# $line"
    done    
    echo "################################################################################"
    echo
}
