#!/bin/sh
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2007 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

set -x

export DISPLAY='none'
export MCR_INHIBIT_CTF_LOCK=1

machine=`uname -m`

if [ "$machine" = "x86_64" ]
then extension=64
else extension=32
fi

thisdir=`dirname $0`

$thisdir/run_birdseye_$extension.sh $*
