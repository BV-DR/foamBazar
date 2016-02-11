#! /bin/bash

if [ $# -lt 3 ]; then
cat << USAGE
Create cellSetDict_box file(s)

Usage: $0 X dX nX [Ymin Ymax] [Zmin Zmax]

    X       start pos.
    dX      increament
    nX      number of slices
    Ymin    default to -20
    Ymax    default to +20
    Zmin    default to -25
    Zmax    default to +25

USAGE
exit;
fi

X=$1; dX=$2; nX=$3; FOUT="cellSetDict_box";
Ymin="-20"; Ymax="20"; Zmin="-25"; Zmax="25"; NAME="box";
if [ $# -ge 4 ]; then Ymin=$4; fi
if [ $# -ge 5 ]; then Ymin=$5; fi
if [ $# -ge 6 ]; then Ymin=$6; fi
if [ $# -ge 7 ]; then Ymin=$7; fi

cat > ".tmp_${FOUT}" << EOF
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.5                                   |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      cellSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

name NAME;

action new;

topoSetSources
(
    boxToCell
    {
        box (X1 Ymin Zmin) (X2 Ymax Zmax);
    }
);

// ************************************************************************* //
EOF

i=0; while [ $i -lt $nX ]; do
X1=$(echo "$X + $dX*$i" | bc -l)
X2=$(echo "$X1 + $dX" | bc -l)
echo "$X1 $X2"
let i=i+1
sed -e "s/X1/$X1/" -e "s/X2/$X2/" -e "s/Ymin/$Ymin/" -e "s/Ymax/$Ymax/" -e "s/Zmin/$Zmin/" -e "s/Zmax/$Zmax/" -e "s/NAME/$NAME$i/" ".tmp_${FOUT}" >  "${FOUT}${i}"
done

rm -f ".tmp_${FOUT}"

