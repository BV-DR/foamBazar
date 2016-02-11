#! /bin/bash

logfile="log.translate0"
if [ ! -f $logfile ]; then
transformPoints -translate "(0 0 -15430.6)" 2>&1 | tee $logfile
fi

logfile="log.snappyHexMesh"
if [ ! -f $logfile ]; then
    snappyHexMesh 2>&1 | tee $logfile
    # time directories in REVERSE order
    timeDir=`ls -d */ | sort -g | tr 'a-d' "*" | tr 'f-z' "*" | grep -v "*" | sed -e "s/\//\n/g" | sort -gr | grep -v "\.\." | sed "/^$/d"`
    if [ -d "shm.d1" ]; then rm -fr shm.d1; fi
    mkdir -p shm.d1
    for i in $time; do if [ "$i" == "0" ]; then continue; fi
        mv $i shm.d1
    done
    for i in $time; do
        yes | cp -r shm.d1/$i/* constant/
    done
fi

logfile="log.translate1"
if [ ! -f $logfile ]; then 
transformPoints -translate "(0 0 15430.6)" 2>&1 | tee $logfile
fi

logfile="log.scale"
if [ ! -f $logfile ]; then
transformPoints -scale "(1e-3 1e-3 1e-3)" 2>&1 | tee $logfile
fi

logfile="log.translate2"
if [ ! -f $logfile ]; then
transformPoints -translate "(0 0 0.82799)" 2>&1 | tee $logfile
fi

