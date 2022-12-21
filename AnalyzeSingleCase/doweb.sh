#!/bin/bash
#
# This script archive the figures and makes
# a web page using the template $fhtml
#
# Qing Li, 160501
#
if [ $# -eq 9 ]; then
	casename=$1
	fhtml=$2
	arcdir=$3
	matfile=$4
	astart=$5
	aend=$6
	utau=$7
	wtsfc=$8
	fcor=$9
else
	echo "Need 9 arguments. Stop!"
	exit 1
fi

#
casedir=$arcdir/profiles/$casename
# create a folder if not exist
mkdir -p $casedir
# move mat file to $casedir
mv $matfile $casedir/
# move *.eps to $casedir
mv *.eps $casedir/
rm *.fig
# convert .eps to .gif
for f in $casedir/*.eps
do
	convert -density 150 $f ${f//eps/gif}
done
# cp html to case directory
cp $fhtml $casedir/index.html
sed -i.bk "s#CASE_NAME#$casename#g" $casedir/index.html
sed -i.bk "s#AVGST#$astart#g" $casedir/index.html
sed -i.bk "s#AVGED#$aend#g" $casedir/index.html
sed -i.bk "s#UTAU#$utau#g" $casedir/index.html
sed -i.bk "s#WTSFC#$wtsfc#g" $casedir/index.html
sed -i.bk "s#FCOR#$fcor#g" $casedir/index.html
rm $casedir/index.html.bk
# archive .eps for later use
mkdir -p $casedir/eps
mv $casedir/*.eps $casedir/eps/
