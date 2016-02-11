#!/bin/bash
for i in `ls -d forces_*`; do
    rm -fr matlab/${i#*_}
    mkdir -p matlab/${i#*_}
    (
    cd ${i};
    for j in `ls -d * | sort -g`; do
    echo ${i}/$j/forces.dat
    sed -e "s/#/%/g" -e "s/(/ /g" -e "s/)/ /g" ${j}/forces.dat >> ../matlab/${i#*_}/forces.dat;
    done
    )
done

exit
for i in `ls -d forceCoeffs_*`; do
    mkdir -p matlab/${i#*_}
    for j in `ls -rt ${i}/*/forceCoeffs.dat | sort -n`; do
    sed -e "s/#/%/g" ${j} >> ./matlab/${i#*_}/forceCoeffs.dat;
    done
done


