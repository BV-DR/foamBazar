#! /bin/bash
if [ $# -lt 1 ]; then
cat << USAGE
Create subsetMesh for every time directory. Cannot handle topology changed yet...
FIXME: make it work with a decomposed case

Usage:
    $0 cellSetDict [nice]
USAGE
exit;
fi

if [ $# -ge 2 ]; then NICE=$2; else NICE=0; fi

cellSetDict=$1
cellSetName=`grep name $cellSetDict | sed -e "s/name //" -e "s/;//"`

# don't need this line for now, patchNames=${*/$1/} # cut away the 1st arg

# prepare for output
targetDir="mksubset.d/`basename $cellSetDict/`"; rm -fr $target
mkdir -p ${targetDir}/{constant,system}
lndir -silent ../../../constant ${targetDir}/constant
lndir -silent ../../../system ${targetDir}/system
(cd $targetDir && ln -s ../../* . > /dev/null 2>&1)
timeDir=`ls -d */ | sort -g | tr 'a-d' "*" | tr 'f-z' "*" | grep -v "*"` # get time directories
for i in $timeDir; do rm -fr ${targetDir}/${i/\//}; done # clear sym.links of timeDir

# create cellSet based on "constant" directory
cp -f $cellSetDict ${targetDir}/system/cellSetDict
nice -n $NICE cellSet -case ${targetDir}

# loop through timeDir, 
for i in $timeDir; do
    echo "processing: $i"
    (cd $targetDir && rm -fr "tmp_${i}") # this is necessary
    mkdir -p "${targetDir}/${i}/polyMesh/"
    lndir ../../../${i} ${targetDir}/${i} > /dev/null
    (cd ${targetDir}/${i}/polyMesh/ && ln -sf ../../constant/polyMesh/sets/ ./sets)
    (cd $targetDir && nice -n $NICE subsetMesh $cellSetName -overwrite)
    (cd $targetDir && mv ${i} "tmp_${i}") # this is necessary
done

for i in $timeDir; do
     (cd $targetDir && mv "tmp_${i}" ${i})
done

