#! /bin/bash

if [ $# -lt 2 ]; then
cat << USAGE

Usage: `basename $0` fieldName patchName(s)

output(s)   patchIntegrate.log
            patchIntegrate.txt (contains "Time Total.Int(x, y, z)" )

USAGE
exit;
fi

fieldName=$1;
patchName=${*//$1/};
logfile="./patchIntegrate_${1}.log"
txtfile="./patchIntegrate_${1}.txt"

if [ ! -f $logfile ]; then
    patchIntegrate -fieldName $fieldName -patchName "($patchName)" > $logfile
else
    echo
    echo "Using existed file: $logfile"
    echo
fi

# grep the total intergral
grep Time $logfile | sed -e "1d" -e "s/Time = //" > .tmp
grep "Total in" $logfile | sed -e "s/    Total int. = (//" -e "s/)//" > .tmp1
if [ `cat .tmp | wc -l` -gt `cat .tmp1 | wc -l` ]; then sed -i -e "1i 0 0 0" .tmp1; fi
paste -d " " .tmp .tmp1 > $txtfile
rm -f .tmp .tmp1

# grep min/max field value(s)
grep "FieldMin" $logfile | sed -e "s/.*=\ //" | awk 'min=="" || $1 < min {min=$1} END{ print min}' > ${1}_MinMax
grep "FieldMax" $logfile | sed -e "s/.*=\ //" | awk 'max=="" || $1 > max {max=$1} END{ print max}' >> ${1}_MinMax





