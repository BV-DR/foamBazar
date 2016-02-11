#! /bin/bash
if [ $# -lt 2 ]; then
cat << EOF

Combine results from "time(s)" dir into one single "dir"

Usage: $0 dir time(s)

EOF
exit;
fi

fdir=$1
shift

for tdir in "$@"; do
if [ "$tdir" = "$fdir" ]; then continue; fi
(
    cd $fdir;
    for iFile in `ls *`; do
        fin=../$tdir/$iFile
        fout=$iFile
        tname=$tdir
        good=`grep "^$tname " $fout | wc -l`
        if [ $good -eq 0 ]; then
            #echo "ignore file (no match): $fin" 
            continue;
        fi
        
        tlast0=`tail -n1 $fout | sed "s/ .*//"`
        tlast1=`tail -n1 $fin | sed "s/ .*//"`
        good=`echo $tlast0'<'$tlast1 | bc -l`
        if [ $good -eq 0 ]; then
            #echo "ignore file (too short): $fin"
            continue;
        fi

        #echo "add file: $tdir/$iFile to $fdir/" 
        newtname=`grep "^$tlast0 " $fin | wc -l`
        if [ $newtname -eq 0 ]; then
            sed -i "/^\b$tname /q" $iFile
        else
            tname=$tlast0
        fi
        sed -n "/^$tname/,$ p" ../$tdir/$iFile >> $iFile
    done
)
done
