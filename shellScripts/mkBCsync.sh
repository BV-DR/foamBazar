#! /bin/bash

if [ $# -lt 1 ]; then
cat << EOF
Sync. BC according to "constant/polyMesh/boundary" 

Usage: $0 BC-file(s)

EOF
exit;
fi

BC=`grep -B2 type constant/polyMesh/boundary | grep -v "\-\-" | sed -e "/{/,/;/d"`
nBC=`grep -B2 type constant/polyMesh/boundary | grep -v "\-\-" | sed -e "/{/,/;/d" | wc -l`

echo
echo "Number of BC(s) found in \"constant\": $nBC" 
echo
for i in $*; do
    echo "Sync.: $i"
    if [ ${i#.*} == "gz" ]; then gunzip $i 2>/dev/null; i=${i%*.gz}; fi
    iBC=`grep -B2 type $i | grep -v "\-\-" | sed -e "/{/,/;/d"`
    ni=`grep -B2 type $i | grep -v "\-\-" | sed -e "/{/,/;/d" | wc -l`
    echo -n "$ni "
    if [ $ni -eq $nBC ]; then echo "...skip"; continue; fi
    echo -n "Remove: "
    for j in $iBC; do
        match=0;
        for jj in $BC; do if [ "$j" == "$jj" ]; then match=1; break; fi; done
        if [ $match -eq 0 ]; then # dosn't match, let's remove it
            echo -n "$j "
            sed -i -e "/${j}/,/}/d" $i  
        fi
    done
    echo
done




