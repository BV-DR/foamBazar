#! /bin/bash
#
# create a subdirectory, links files, change BC, decompose, restore BC in processor directories, exit
#
# the processor* directories are the final product.
#
if [ "$1" == "-help" ]; then
cat << EOF

decomposePar is sick when dealing with nonuniform pointField. The solution is to define the boundaries to be uniform and decompose the case. After the case has been decomposed, the correct boundary condition can be applied in each of the processor-directory.

Usage: `basename $0` boundary-file target-file(s)

EOF
exit
fi
if [ ! -d "0" ]; then echo; echo "0 time directory does not exist... terminate."; echo; exit 1; fi
if [ ! -d "constant" ]; then echo; echo "constant directory does not exist... terminate."; echo; exit 1; fi
if [ ! -d "system" ]; then echo; echo "system directory does not exist... terminate."; echo; exit 1; fi

# prepare the subdirectory
subdir=`basename $0`; subdir=${subdir%*.sh}.d; rm -fr $subdir
mkdir -p $subdir
cd $subdir
ln -s ../{constant,system} .
cp -r ../0 .
for i in `ls 0/*.gz 2> /dev/null `; do yes | gunzip $i; done
cp -r 0 0-orig
sed -i -e "s/\$FOAM_CASE\//MOIN_FOAM_CASE_MOIN/g" 0-orig/*

# change all patch(es) to "type slip;"
patch=`grep -B2 type constant/polyMesh/boundary | grep -v "\-\-" | sed -n 'p;N;N'`
for jj in `ls 0/*`; do
    sed -i -e "/\<type\>/,/}/d" $jj
    for i in $patch; do
        sed -i -e "/\<$i\>/,+1 s/{/{\n        type slip;\n    }/"  $jj
    done
done

#undopatch=`grep -B2 type constant/polyMesh/boundary | grep -v "\-\-" | grep -v "patch" | grep -v "wall" | grep -B2 "type" | grep -v "\-\-" | sed -n 'p;N;N'`

# restore those with type diff.from "wall" or "patch"
undopatch=`grep -B2 type constant/polyMesh/boundary | grep -v "\-\-" | grep -v "patch" | grep -v "wall" | grep -B2 "type" | grep -v "\-\-" | sed -n 'p;N;N'`
for jj in `ls 0/*`; do
    for i in $undopatch; do
        tmp="0-orig/`basename $jj`"; tmp=`sed -n -e "/\<$i\>/,/}/ p" $tmp | sed -e "1s/\<$i\>//" -e "s/{//" -e "s/}//"`
        tmp=`echo $tmp`
        sed -i -e "/\<$i\>/,/}/ s/type slip;/$tmp/"  $jj
    done
done

decomposePar
 if [ $? -ne 0 ]; then exit; fi

echo
echo -n "Restore BC ... "

for jj in `ls processor*/0/*.gz 2> /dev/null`; do yes | gunzip $jj; done

patch=`grep -B2 type constant/polyMesh/boundary | grep -v "\-\-" | egrep -B2 "wall|patch" | grep -v "\-\-" | sed -n 'p;N;N'`
for jj in `ls processor*/0/*`; do
    for i in $patch; do
        echo "proc.$jj"
        tmp="0-orig/`basename $jj`"; tmp=`sed -n -e "/\<$i\>/,/}/ p" $tmp | sed -e "s/\<$i\>//" -e "s/{//" -e "s/}//" -e "s/\/\/.*//g"`
        tmp=`echo $tmp`
        sed -i -e "/\<$i\>/,/}/ s/type.*slip;/$tmp/"  $jj
    done
done

sed -i -e "s/MOIN_FOAM_CASE_MOIN/\$FOAM_CASE\//g" 0-orig/* processor*/0/*

echo "done"
echo
echo "output in \"$subdir\""
echo

