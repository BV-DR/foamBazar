#!/bin/bash

if [ $# -lt 1 ]; then
cat << EOF
Extract Courant number from log file and make a timeserie plot using gnuplot.

Usage: $0 run.log

EOF
exit;
fi



grep "^Time =\|Courant" $1 | grep -v "Mesh" | awk 'BEGIN{RS="max: "} {print$4,$1}' | sed -e "1,3d" | grep -v ":" > Courant.dat
#sed -i -e "1,3d" Courant.dat

j=0; for i in `tail -n1 Courant.dat`; do let j=j+1; done;
if [ $j -eq 1 ]; then sed -i '$d' Courant.dat; fi

if [ $# -lt 2 ]; then
    gnuplot -persist -e 'plot "Courant.dat" with lp;'
    exit;
fi

if [ $2 == "-" ]; then
    cat Courant.dat
    exit;
fi

if [ $2 == "r" ]; then 
if [ $# -lt 3 ]; then PAUSE=5; else PAUSE=$3; fi
if [ $# -ge 4 ]; then 
    j=0;
    for i in $*; do let j=j+1;
        if [ $j -ge 4 ]; then
            plotopt="$plotopt $i"
        fi    
    done
    echo $plotopt 
fi
tmpfile=/tmp/gnuplot_$RANDOM.gpl
cat > $tmpfile << EOF
plot "< $0 $1 -" $plotopt 
pause $PAUSE
replot
reread
EOF
gnuplot -persist $tmpfile 
exit
fi

