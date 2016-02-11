#! /bin/bash
if [ $# -lt 1 ]; then echo "Usage: $0 patchname(s)"; exit; fi


FID="controlDict_mkforces"; if [ -f $FID ]; then rm -f $FID; fi
cp -f system/controlDict $FID 2> /dev/null

echo >> $FID
echo -e "functions\n(" >> $FID


for i in $*; do
cat >> $FID << EOF
    forces_$i
    {
        type    forces;
        functionObjectLibs ("libforces.so");
        patches ($i);
        rhoInf  1000.0;
        CofR    (0 0 0);
    }
EOF
done

echo ");" >> $FID


