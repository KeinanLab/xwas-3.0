#! /bin/bash
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2007 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

# If there are two executables in the same directory foo, and foo.64, then symlink
# this script into the same directory with the name foo.sh.  It will detect whether to
# run 32 or 64-bit, and will launch the appropriate program with all args.

machine=`uname -m`

if [[ "$machine" = x86_64 ]]
then extension=.64
else extension=
fi


thisdir=`dirname $0`
progname=`basename $0 .sh`

$thisdir/$progname$extension $*
