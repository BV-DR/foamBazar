#!/bin/bash

if [ $# -lt 1 ]; then
cat << EOF
Extract (min/max) gamma from log file and make a timeserie plot using gnuplot.

Usage: $0 run.log

EOF
exit;
fi


grep "^Time =" $1  | sed -e "s/Time = //" > timegamma.dat
grep "Max(gamma)" $1 | sed -e "s/.*Max(gamma) = //" > maxgamma.dat
grep "Max(gamma)" $1 | sed -e "s/.*Min(gamma) = //" -e "s/ Max.*//" > mingamma.dat
paste timegamma.dat mingamma.dat maxgamma.dat > gamma.dat
rm -f timegamma.dat mingamma.dat maxgamma.dat
gnuplot -persist -e 'plot "gamma.dat" using 1:2 with l, "gamma.dat" using 1:3 with l;'
