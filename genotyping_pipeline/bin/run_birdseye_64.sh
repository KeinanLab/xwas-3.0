#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
# NOTE: THIS IS A HACKED VERSION OF THE FILE PRODUCED BY MATLAB COMPILER
# This script contains the following changes from the one produced by Matlab compiler:
# * Assumes that MCR (Matlab compiled runtime) resides in a subdirectory
# of the directory containing this script.
# Does not require that the birdseye executable be in the PWD.

DIRNAME=`dirname $0`
DIRNAME=`${DIRNAME}/getAbsPath ${DIRNAME}`

MCRROOT=$DIRNAME/MCR75_glnxa64/v77
echo "------------------------------------------"
  echo Setting up environment variables
  echo ---
  MWE_ARCH="glnxa64" ;
  if [ "$MWE_ARCH" = "sol64" ] ; then
	LD_LIBRARY_PATH=.:/usr/lib/lwp:${MCRROOT}/runtime/glnxa64 ; 
  else
  	LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  fi
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  if [ "$MWE_ARCH" = "mac" -o "$MWE_ARCH" = "maci" ]; then
	DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/System/Library/Frameworks/JavaVM.framework/JavaVM:/System/Library/Frameworks/JavaVM.framework/Libraries;
  else
	MCRJREVER=`cat ${MCRROOT}/sys/java/jre/glnxa64/jre.cfg` ; 
	echo Found MCR Java JRE version: $MCRJREVER ; 
	MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre${MCRJREVER}/lib/amd64 ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
  fi
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export LD_LIBRARY_PATH;
  export XAPPLRESDIR;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
$DIRNAME/birdseye_64 $*
exit
