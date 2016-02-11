#! /bin/bash
if [ $# -lt 1 ]; then echo "mkcellZones.sh: create \"cellZones\" for viewing in                                                        paraview"; exit 1; fi
if [ -f "./cellZones.gz" ] | [ -f "cellZones" ]; then
        rm -f cellZones*
fi

cat > .tmp << EOF
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.5.x                                 |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       regIOobject;
    location    "";
    object      cellZones;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

EOF
sed -e 's/\\/\\\\/g' .tmp

echo -e "$#\n(" >> .tmp
for i in $*; do
        bn=`basename $i`
        echo -e "$bn\n{\n\ttype cellZone;\ncellLabels\tList<label>" >> .tmp
        sed -e "1,18d" $i | sed '$d' | sed '$d' >> .tmp
        echo -e ";\n}" >> .tmp
done
echo -e ")\n\n// ***************************************************************                                                       ********** //" >> .tmp
mv -f  .tmp cellZones

