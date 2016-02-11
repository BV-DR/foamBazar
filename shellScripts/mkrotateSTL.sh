#! /bin/bash

# return VERTEX
# usage: matrix_vector_add x y z
# info: subtract a vector(x,y,z) to a matrix of matrix(x,y,z, x,y,z, ...)
function move_to_cog {
j=1
for i in $VERTEX; do
    if [ $j -eq 1 ]; then tmp = "$tmp `echo $i - $1 | bc -l`"; j=2;
    else if [ $j -eq 2 ]; then tmp = "$tmp `echo $i - $2 | bc -l`"; j=3;
    else tmp = "$tmp `echo $i - $3 | bc -l`"; j=1;
    fi fi
done
VERTEX=$tmp
}
function move_back_from_cog {
j=1
for i in $VERTEX; do
    if [ $j -eq 1 ]; then tmp = "$tmp `echo $i + $1 | bc -l`"; j=2;
    else if [ $j -eq 2 ]; then tmp = "$tmp `echo $i + $2 | bc -l`"; j=3;
    else tmp = "$tmp `echo $i + $3 | bc -l`"; j=1;
    fi fi
done
VERTEX=$tmp
}

# return SOLID, FACE, VERTEX
# usage: read_STL filename
# info: - 
function read_STL {
SOLID=`grep solid $1 | grep -v end | sed -e "s/solid //"`
FACE=`grep "facet normal" $1 | sed -e "s/facet normal //"`
VERTEX=`grep vertex $1 | sed -e "s/    vertex //"`
}
function write_STL {
fout=$1; rm -f $1;
echo "solid $SOLID" >> $fout
while [ -n "$FACE" ]; do
    echo -n "facet normal " >> $fout
    jj=1
    for j in $FACE; do
        tmp="$tmp $j"
        if [ $jj -eq 3 ]; then break; else let jj=jj+1; fi
    done
    echo $tmp >> $fout
    FACE=${FACE//$tmp/}
    tmp=""
    echo "endfacet" >> $fout
done
echo "endsolid $SOLID" >> $fout
}

# rotate FACE, VERTEX
function rotate_STL {
Te[1]=`echo "c($3)*c($2)" | bc -l`;
Te[2]=`echo "c($3)*s($2)*s($1)-s($3)*c($1)" | bc -l`;
Te[3]=`echo "c($3)*s($2)*c($1)+s($3)*s($1)" | bc -l`;
Te[4]=`echo "s($3)*c($2)" | bc -l`;
Te[5]=`echo "s($3)*s($2)*s($1)+c($3)*c($1)" | bc -l`;
Te[6]=`echo "s($3)*s($2)*c($1)-c($3)*s($1)" | bc -l`;
Te[7]=`echo "-s($2)" | bc -l`;
Te[8]=`echo "c($2)*s($1)" | bc -l`;
Te[9]=`echo "c($2)*c($1)" | bc -l`;
#[cosd(zdeg)*cosd(ydeg) cosd(zdeg)*sind(ydeg)*sind(xdeg)-sind(zdeg)*cosd(xdeg) cosd(zdeg)*sind(ydeg)*cosd(xdeg)+sind(zdeg)*sind(xdeg)
#      sind(zdeg)*cosd(ydeg) sind(zdeg)*sind(ydeg)*sind(xdeg)+cosd(zdeg)*cosd(xdeg) sind(zdeg)*sind(ydeg)*cosd(xdeg)-cosd(zdeg)*sind(xdeg)
#                -sind(ydeg) cosd(ydeg)*sind(xdeg) cosd(ydeg)*cosd(xdeg)];
j=1
for i in $VERTEX; do
    tmp1=`echo "$tmp1+$Te[$j]*$i" | bc -l`
    tmp2=`echo "$tmp2+$Te[$j+3]*$i" | bc -l`
    tmp3=`echo "$tmp3+$Te[$j+6]*$i" | bc -l`
    if [ $j -eq 3 ]; then 
        j=1;
        tmp="$tmp $tmp1 $tmp2 $tmp3"
        tmp1=0; tmp2=0; tmp3=0;
    else
        let j=j+1
    fi
done
VERTEX=$tmp
}

read_STL $1
#move_to_cog $cogx $cogy $cogz
#rotate_STL $xdeg $ydeg $zdeg 
#move_back_from_cog $cogx $cogy $cogz
#translate_MWL $MWL
write_STL ${1}.stl


 
