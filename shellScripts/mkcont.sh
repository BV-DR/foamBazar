#! /bin/bash
if [ $# -lt 1 ]; then
cat << EOF

Restart a previously paused openfoam JOB (parallel or serial run)

Usage: $0 run.log
       
       a log file is required to lookup the pid(s)

EOF
exit;
fi

logfile=$1
nProcs=`head -n25 $logfile | grep nProcs | sed -e "s/nProcs.*: //"`
PID=`head -n25 $logfile | grep PID | sed -e "s/PID.*: //"`
tmp=`head -n$[ $nProcs + 25 ] $logfile | grep -A$[$nProcs + 2] Slaves | grep "\." | sed -e "s/\"//g;s/.*\.//" `
PID="$PID `echo $tmp`"

    for i in $PID; do
        echo "Restart pid: $i"
        kill -s CONT $i
    done

