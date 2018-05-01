#!/bin/sh

# Liftover coordinates and creates a new map file with new coordinates. 
fname=$1
wname=$2
build=$3

if [ -f $wname ]; then
    rm $wname
fi

if [ -f tomap ]; then
    rm tomap
fi

echo "Reading $fname into tomap file..."

while read line
do 
	chr=$(echo ${line} | awk '{print $1}' )
	if [ $chr = '23' ] || [ $chr = 'X' ]
	then
		chr="23"
		echo ${line} | awk '{print "chrX",$4,$4+1,$2}'  >> tomap
	else
		echo ${line} | awk '{print "chr"'$chr',$4,$4+1,$2}'  >> tomap
	fi
done <"$fname"

echo "Calling liftOver..."
./liftOver tomap hg${build}ToHg19.over.chain mapsuc mapfail

if [ -f mapkey ]; then
    rm mapkey
fi

sz=$(du -sk mapfail | awk '{print $1}')
if [ ${sz} -gt 0 ]
then
	cat mapfail | grep -v new | awk '{print $4,0}' >> mapkey
fi
awk '{print $4,$2}' ./mapsuc >> mapkey

awk 'FNR==NR {f1[$1]=$2; next} ($2 in f1) {print $1,"\t"$2"\t0\t"f1[$2]}' ./mapkey ./$fname > $wname

rm tomap mapsuc mapfail mapkey
