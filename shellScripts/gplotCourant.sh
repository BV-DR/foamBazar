#!/bin/bash

if [ $# -lt 1 ]; then
cat << EOF
Extract Co. number from log file and make a timeserie plot using gnuplot.

Usage: $0 run.log #remove_n_lines

EOF
exit;
fi

logfile=$1
title="`basename $(pwd)` : $logfile"
if [ $# -ge 2 ]; then removeLastN=$2; else removeLastN=0; fi

cmd1="<cat $logfile | grep -w -E \"^Time =|^Courant Number\" | grep -B1 \"^Time\" | sed \"s/--//\" | sed \"/^\$/d\" | sed \"/^Time =/N;s/\n/ /g\" | sed \"s/Time = //;s/Courant Number mean: / /;s/max: / /;s/velocity magnitude:/ /\" | sed \"1d;\$d\" | head -n -$removeLastN "

cmd2="<cat $logfile | grep -w -E \"^Time =|^Interface Courant Number\" | grep -B1 \"^Time\" | sed \"s/--//\" | sed \"/^\$/d\" | sed \"N;/^Interface Courant Number/s/\n/ /\" | sed \"s/Time = //;s/^Interface Courant Number mean: //;s/max: / /\" | sed \"1d;\$d\" | head -n -$removeLastN "

gnuplot -persist -raise -e "set title '$title'; plot \
'$cmd1' u 1:2 w l t \"Co. (mean)\" \
, '$cmd1' u 1:3 w l t \"Co. (max)\" \
, '$cmd2' u 3:1 w l t \"Interface Co. (mean)\" \
, '$cmd2' u 3:2 w l t \"Interface Co. (max)\" \
"


