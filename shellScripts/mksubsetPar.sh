#! /bin/bash
if [ $# -lt 1 ]; then
cat << USAGE
Create subsetMesh for every time directory. Cannot handle topology changed yet...
NOTE: For a decomposed case only

Usage:
    $0 cellSetDict [nice]
USAGE
exit;
fi
if [ ! -d processor0 ]; then echo; echo "${0}: processor* directory(s) not found" ; echo; exit; fi
if [ $# -ge 2 ]; then NICE=$2; else NICE=0; fi

cellSetDict=$1
cellSetName=`grep name $cellSetDict | sed -e "s/name //" -e "s/;//"`
np=`ls -d processor* | wc -l`

# don't need this line for now, patchNames=${*/$1/} # cut away the 1st arg

# prepare for output
targetDir="mksubset.d/`basename $cellSetDict/`"; rm -fr $targetDir
mkdir -p ${targetDir}/{constant,system}
lndir -silent ../../../constant ${targetDir}/constant
lndir -silent ../../../system ${targetDir}/system
j=0; while [ $j -lt $np ]; do
    mkdir -p ${targetDir}/processor${j}/constant
    lndir -silent ../../../../processor${j}/constant ${targetDir}/processor${j}/constant
    let j=j+1
done
(cd $targetDir && ln -s ../../* . > /dev/null 2>&1)
timeDir=`ls -d processor0/*/ | sort -g | sed -e "s/processor0\///" | tr 'a-d' "*" | tr 'f-z' "*" | grep -v "*"` # get time directories
for i in $timeDir; do rm -fr ${targetDir}/${i/\//}; done # clear sym.links of timeDir

# create cellSet based on "constant" directory
cp -Lf $cellSetDict ${targetDir}/system/cellSetDict
nice -n $NICE mpirun -np $np cellSet -case ${targetDir} -parallel
echo

# loop through timeDir, 
for i in $timeDir; do
    echo "processing: $i"
    (cd $targetDir && rm -fr "processor*/tmp_${i}") # this is necessary
    j=0; while [ $j -lt $np ]; do
        mkdir -p "${targetDir}/processor${j}/${i}/polyMesh/"
        lndir ../../../../processor${j}/${i} ${targetDir}/processor${j}/${i} > /dev/null
        (cd ${targetDir}/processor${j}/${i}/polyMesh/ && ln -sf ../../constant/polyMesh/sets/ ./sets)
        let j=j+1
    done
    (cd $targetDir && nice -n $NICE mpirun -np $np subsetMesh $cellSetName -overwrite -parallel)
    j=0; while [ $j -lt $np ]; do
        (cd $targetDir/processor${j} && mv ${i} "tmp_${i}") # this is necessary
        let j=j+1
    done
done

for i in $timeDir; do
    j=0; while [ $j -lt $np ]; do
        (cd ${targetDir}/processor${j} && mv "tmp_${i}" ${i})
        let j=j+1
    done
done

