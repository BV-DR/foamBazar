#! /bin/bash
if [ $# -lt 1 ]; then
cat << EOF
Extract ExecutionTime of each timestep from log file and make a timeserie plot using gnuplot.

Usage: $0 run.log N(default,N=1)

EOF
exit;
fi

if [ $# -ge 2 ]; then N=$2; else N=1; fi

grep "^Time =" $1  | sed -e "s/Time = //" | sed -n -e "0~${N}p;" > timeCPUtime.dat
grep ExecutionTime $1 | sed -n -e "0~${N}p;" | awk '{y=x; x=$3}{print x-y}' > diffCPUtime.dat
paste timeCPUtime.dat diffCPUtime.dat > CPUtime.dat
i=`cat timeCPUtime.dat | wc -l`; j=`cat diffCPUtime.dat | wc -l`;
if [ $i -ne $j ]; then sed -i '$d' CPUtime.dat; fi
#rm -f timeCPUtime.dat diffCPUtime.dat
gnuplot -persist -e 'plot "CPUtime.dat" with lp'

