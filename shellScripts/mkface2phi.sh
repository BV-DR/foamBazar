#! /bin/bash

fid='0/faceOrthogonality';
fbc='constant/polyMesh/boundary';

if [ -n "grep scalarList $fid" ]; then echo "$fid: not a scalaList() ...exit"; exit; fi

cp -f $fid /tmp/`basename $fid`
cp -f $fbc /tmp/`basename $fbc`

fid0=$fid
fid=/tmp/`basename $fid`
fbc=/tmp/`basename $fbc`

nFaces=`sed -n "19p" $fid`

BC=`grep -B1 "^    {" $fbc | grep -v -E "\-\-|{"`
nBC=`grep -B1 "^    {" $fbc | grep -v -E "\-\-|{"|wc -l`
idx0=20
j=-1;
for i in $BC; do
    startf=`sed -n "/$i/,/^    }/p" $fbc|grep startFace|sed "s/;//;s/^.* //"`
    n=`sed -n "/$i/,/^    }/p" $fbc|grep nFaces|sed "s/;//;s/^.* //"`
    let j=j+1; idx[$j]=`expr $idx0 + $startf`
    echo $startf $n ${idx[$j]} `expr $n + ${idx[$j]}`
done
    let j=j+1; idx[$j]=`expr $idx0 + $startf + $n`

i=0; sed -i "${idx[$i]} s/\$/); boundaryField { /" $fid

for ii in $BC; do
    sed -i "${idx[$i]} s/\$/$ii {type calculated;value uniform 0;} /" $fid
done

i=0; sed -i "${idx[$i]} s/\$/}\n/" $fid

sed -i "`expr 2 + ${idx[$i]}`,\$ d" $fid

sed -i "
    12 s/scalarList;/surfaceScalarField;/
    18 i dimensions [0 0 0 0 0 0 0];
    19 i internalField nonuniform List<scalar>
    19 s/^.*/`expr ${idx[0]} - $idx0`/
" $fid

mv $fid $fid0

