#! /bin/bash
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2007 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

set -e

# If --exeDir has been passed on the command line, it must be detected so that
# R_LIBS can be set correctly, in addition to passing --exeDir along to Birdsuite.jar.
# If --exeDir has not been passed on the command line, then by default the
# executables are found in the same directory as birdsuite.sh itself.

exedir=`dirname $0`
exedir_implicit=1

# This mess of punctuation iterates downward, because
# the command-line args go into BASH_ARGV backward.
# If --exeDir is the last argument on the command line, it is
# not detected here, but will cause Birdsuite.jar to fail.
for (( i=((${#BASH_ARGV[@]}-1)) ; i > 0 ; --i ))
do  if [[ ${BASH_ARGV[$i]} == --exeDir ]]
    then exedir=${BASH_ARGV[(($i-1))]}
        exedir_implicit=''
    fi
done

set -x

unset DISPLAY
export MCR_INHIBIT_CTF_LOCK=1
#experimental change for when python is packaged.
#PYTHONPATH=$PYTHONPATH:$exedir
#export PYTHONPATH
exedir=`${exedir}/getAbsPath ${exedir}`

# So that R can pick up broadgap packages
export R_LIBS=$exedir:$R_LIBS

if [ -n "$exedir_implicit" ]
then exedir_args="--exeDir=$exedir"
fi

java -jar $exedir/Birdsuite.jar $exedir_args $* --exeDir=$exedir --metadataDir=$exedir/../metadata
