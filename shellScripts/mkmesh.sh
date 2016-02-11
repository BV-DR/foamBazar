#! /bin/bash

if [ $# -lt 1 ]; then
cat << EOF

Create openfoam mesh from gmsh geo-file.
	
	usage: `basename $0` file.geo [gmshopt(s)]

EOF
exit;
fi

geofile=$1;

j=0; for i in $*; do if [ $j -eq 0 ]; then let j=j+1; continue; fi; gmshopts="$gmshopts $i";  done

rm -fr constant/polyMesh
(
    gmsh $gmshopts -3 -format msh2 $geofile 
    gmshToFoam ${geofile%*.geo}.msh
    mkSetGmshBC.sh constant/polyMesh/boundary
    rm0FacePatch.sh constant/polyMesh/boundary i
    renumberMesh -overwrite -constant 
    checkMesh
) 2>&1 | tee log.mkmesh



