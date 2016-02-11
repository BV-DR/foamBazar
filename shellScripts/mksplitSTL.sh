#!/bin/bash
if [ $# -lt 1 ]; then
cat << EOF
split a stl-file with multiple solids into multiple one solid stl-files

Usage: $0 stlfile

EOF
exit;
fi

stlfile=$1

for i in `grep -w solid $stlfile | sed -e "s/solid //" -e "s/\r//"`; do
    sed -n -e "/\<solid ${i}\>/,/\<endsolid ${i}\>/p" $stlfile > ${i}.stl
    echo "extract: $i.stl"
done


