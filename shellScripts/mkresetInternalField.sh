#! /bin/bash

if [ $# -lt 1 ]; then
cat << EOF
reset internalField 

Usage: $0 file(s)

EOF
exit;
fi

while [ $# -gt 0 ]; do
    fid=$1; ext=${fid##*.}
    if [ $ext = "gz" ]; then gunzip $fid 2>&1 > /dev/null; fid=${fid%.gz}; fi
    lines=`grep -A1 internalField $fid | sed "1d" | sed "s/^.* //"`
    let lines=lines-1
    if [ $lines -lt 0 ]; then  shift; continue; fi
    sed -i -e "/internalField/,/;/d" $fid
    class=`sed -n -e "/FoamFile/,/}/p" $fid | grep class | sed -e "s/class//" -e "s/.* //" -e "s/;.*//"`;
    case $class in
    "volScalarField")
        type="uniform 0;"
        ;;
    "volVectorField")
        type="uniform (0 0 0);"
        ;;
    esac
    sed -i -e "20 i internalField   $type" $fid
    if [ $ext = "gz" ]; then gzip $fid 2>&1 > /dev/null; fi
    echo -e "reset internalField: $fid"
shift
done

