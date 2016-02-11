#!/bin/bash
#
# just delete every second lines in force*.dat files to reduce the data
#

SRCDIR="matlab"; if [ -d 'matlab-sed' ]; then SRCDIR="matlab-sed"; fi
DESTDIR="matlab-sed"
TMPDIR="mkreduce.tmp"
for ID in `ls -d ${SRCDIR}/*`; do ID=`basename $ID`;
    mkdir -p ${TMPDIR}/${ID}
    sed "n;d" ${SRCDIR}/${ID}/forces.dat > ${TMPDIR}/${ID}/forces.dat;
    sed "n;d" ${SRCDIR}/${ID}/forceCoeffs.dat > ${TMPDIR}/${ID}/forceCoeffs.dat;
done

mkdir -p ${DESTDIR}
cp -rf ${TMPDIR}/* ${DESTDIR}
rm -fr ${TMPDIR}
