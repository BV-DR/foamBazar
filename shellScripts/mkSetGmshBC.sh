#! /bin/bash

if [ $# -lt 1 ]; then
cat << EOF
Change BC type according to Gmsh patch(s)

Usage: `basename $0` boundaryfile
     
     e.g:  body_wall --> patchname is "body" and the type is "wall"
           top_symmetryPlane --> patchname is "top" and the type is "symmetryPlane"

EOF
exit;
fi

filename=$1

patchToChange=`grep -B2 type $filename | grep -v "\-\-" | egrep -B2 "wall|patch" | grep -v "\-\-" | sed -n 'p;N;N' | grep "_"`

for patch in $patchToChange; do
    name=${patch#*_};
    echo "Changing: $patch  -->  $name"
    sed -i -e "/\b${patch}/,/}/s/patch;/${name};/" -e "/\b${patch}/,/}/s/_${name}//" $filename
done

##
###########################################
# the following is a just print out of the current BC  -  not necessary
echo
patches=`grep -B2 type $filename | grep -v "\-\-"  | grep -v "\-\-" | sed -n 'p;N;N'`
for patch in $patches; do
    echo `sed -n -e "/\b${patch}/,/}/p" $filename`
done
echo
