#! /bin/bash

if [ $# -lt 1 ]; then
cat << USAGE
Create mapFields for each time directory. 

Usage:
    $0  <source dir>  [-parallelSource] 

USAGE
exit;
fi

# prepare for output
timeDir=`ls -d */ | sort -g | tr 'a-d' "*" | tr 'f-z' "*" | grep -v "*" | sed -e "s/\//\n/g" | sort -g | grep -v "\.\." | sed "/^$/d"` # get time directories

for i in $timeDir; do 
    if [ "$i" == "0" ]; then continue; fi
    mv $i wait_$i; 
done

mkdir -p log.mapFields

for i in $timeDir; do
    if [ "$i" == "0" ]; then continue; fi
    mv wait_$i $i
    mapFields $* -sourceTime $i 2>&1 | tee log.mapFields/mapFields.log_${i}
done

