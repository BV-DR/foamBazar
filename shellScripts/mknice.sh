#! /bin/bash
if [ $# -lt 1 ]; then
cat << EOF

Change nice number of a running openfoam case. [nice] is default to 10.

Usage: $0 run.log [nice]

       a log file is required to lookup the pid(s)

EOF
exit;
fi
if [ $# -lt 2 ]; then NICE=10; else NICE=$2; fi

#logfile=$1
#nProcs=`head -n150 $logfile | grep "^nProcs" | sed -e "s/nProcs.*: //"`; echo $nProcs
#PID=`head -n150 $logfile | grep PID | sed -e "s/PID.*: //"`
#tmp=`head -n$[ $nProcs + 25 ] $logfile | grep -A$[$nProcs + 2] Slaves | grep "\." | sed -e "s/\"//g;s/.*\.//"`
#PID="$PID `echo $tmp`"

logfile=$1
nProcs=`grep "^nProcs   :" $logfile | head -n1 | sed -e "s/nProcs.*: //"`;
PID=`grep "^PID      :" $logfile | head -n1 | sed -e "s/PID.*: //"`;
SLAVES=`sed -n "/Slaves/,/^)\$/p" $logfile | grep "\." | sed -e "s/\"//g;s/.*\.//"`
PID="$PID `echo $SLAVES`"

    for i in $PID; do
        sudo renice $NICE $i;
    done

