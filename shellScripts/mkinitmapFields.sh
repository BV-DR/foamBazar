#! /bin/bash
if [ $# -lt 2 ]; then
cat << EOF
Init flow field using mapFields. Start with a very small timestep and map/re-map back into 0-dir

Usage: $0 BASEDIR TIMEDIR COUNT

    BASEDIR     system/mapFieldsDict must exist and system/controlDict shall be prepared with 

                endTime nextWrite;
                writeControl timeStep; 
                writeTime 1; 

    TIMEDIR     Start from a small time step
    COUNT       default is 10, map/remap for 10 time(s)

EOF
exit;
fi

BASEDIR=$1; #'initmapFields'
TIMEDIR1=$2; # '0.0001';
if [ $# -ge 3 ]; then COUNTDOWN=$3; else COUNTDOWN=10; fi

APP=`grep application $BASEDIR/system/controlDict | sed -e "s/application//" -e "s/;//"`; 
echo "application $APP;"

if [ -d "$BASEDIR/processor0" ]; then NP=`ls -d processor* | wc -l`; else NP=0; fi

COUNT=0;
COUNT=`ls -d tmp*_$BASEDIR 2>&1 | grep -v "No such file or directory" | wc -l`;
let COUNT=COUNT-1

while [ $COUNT -lt $COUNTDOWN ]; do let COUNT=COUNT+1
    CURDIR="tmp${COUNT}_$BASEDIR"
    rm -fr $CURDIR > /dev/null; mkdir -p $CURDIR; lndir ../$BASEDIR $CURDIR
(   cd $CURDIR;
    if [ $NP -ne 0 ]; then
        echo; echo -n "initmapFields.sh: $COUNT ..."
        if [ 0 -ne $COUNT ]; then 
            mapFields ../tmp$[$COUNT - 1]_$BASEDIR -parallelSource -parallelTarget -sourceTime $TIMEDIR1 > log.mapFields
            for i in `ls -d processor*/0/pointMo*unmap*`; do mv $i ${i%%.*}; done
            echo -n " done mapping ..."
        fi
        mpirun -np $NP $APP -parallel > run.log;
        mkdir $TIMEDIR1;
        echo "done running"
    else
        echo "`basename $0`: run in serial? not yet implemented .. exit"
        exit
    fi
)
done

let COUNT=COUNT+1
CURDIR="tmp${COUNT}_$BASEDIR"
rm -fr $CURDIR > /dev/null; mkdir -p $CURDIR; lndir ../$BASEDIR $CURDIR
(
    cd $CURDIR;
    if [ $NP -ne 0 ]; then
        echo; echo "initmapFields.sh: $COUNT ... finishing ..."
        mapFields ../tmp$[$COUNT - 1]_$BASEDIR -parallelSource -parallelTarget -sourceTime $TIMEDIR1 > log.mapFields
        for i in `ls -d processor*/0/pointMo*`; do mv $i ${i%%.*}; done
    fi
)

