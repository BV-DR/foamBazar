#! /bin/bash

if [ $# -lt 2 ]; then
cat << USAGE

NOTE: For decomposed cases only

Usage: `basename $0` fieldName patchName(s)

output(s)   patchIntegrate.log
            patchIntegrate.txt (contains "Time Total.Int(x, y, z)" )

USAGE
exit;
fi
if [ ! -d processor0 ]; then echo; echo "${0}: processor* directory(s) not found" ; echo; exit; fi

np=`ls -d processor* | wc -l`
fieldName=$1;
patchName=${*//$1/};

if [ ! -f "./patchIntegrate.log" ]; then
    mpirun -np $np patchIntegrate -parallel -fieldName $fieldName -patchName "($patchName)" > ./patchIntegrate.log
else
    echo
    echo "Using existed file: patchIntegrate.log"
    echo
fi

grep Time patchIntegrate.log | sed -e "1,2d" -e "s/Time = //" > .tmp
grep "Total in" patchIntegrate.log | sed -e "s/    Total int. = (//" -e "s/)//" > .tmp1
paste -d " " .tmp .tmp1 > patchIntegrate.txt
rm -f .tmp .tmp1

