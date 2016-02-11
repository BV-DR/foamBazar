#! /bin/bash

fid="system/funkySetFieldsDict.decomposePar"
fin="system/decomposeParDict"

function usage
{
cat << EOF

Create $fid for a manual decomposition

Usage: $0 [optional filename]

    if a filename is passed as the first argument we jump
    to convert_to_labelList() and exit

    content of file system/decomposeParDict:

        manualCoeffs
        {
            proc.number xmin xmax ymin ymax zmin zmax
        }

        e.g: for 4 cores

        manualCoeffs
        {
            dataFile        "funkyCellDist";
            proc.0 0 3 0 1 -1 1;
            proc.1 3 4 0 1 -1 1;
            proc.2 4 7 0 1 -1 1;
            proc.3 7 10 0 1 -1 1;
        }

EOF
exit
}

function create_dict
{
cat > $fid << EOF
FoamFile
{
    version         2.0;
    format          ascii;
    root            "";
    case            "";
    instance        "system";
    local           "";

    class           dictionary;
    object          funkySetFieldsDict;
}

expressions
(
    reset_everything
    {
        field funkyCellDecomp;
        create true;
        expression "0";
        dimension [0 0 0 0 0 0 0];
    }
);
EOF
return
}

function add_proc
{
if [ $# -lt 7 ]; then return; fi
echo $*
sed -i '$d' $fid
cat >> $fid << EOF
    proc.$1
    {
        field funkyCellDecomp; keepPatches true;
        expression "$1";
        condition "pos().x>=$2 && pos().x<$3 && pos().y>=$4 && pos().y<$5 && pos().z>=$6 && pos().z<$7";
    }
);
EOF
return
}

function execute_funky
{
    CMD='funkySetFields -latestTime -dictExt decomposePar'
    eval $CMD
    return
}

function print_count { printf "%3d %9d %9d %9d\n" $2 $1 $3 $4; }
function convert_to_labelList
{
fid=$1
dict="system/decomposeParDict"
fout="constant/`sed -n -e '/manualCoeffs/,/^}/p' $dict | grep dataFile | sed -e 's/^[ \tdataFile]*//' -e 's/\\\"//g' -e 's/;.*//'`"

gunzip $fid 2>&1 | grep -v gzip
fid=${fid%.gz}

sed -e "/FoamFile/,/}/s/volScalarField/labelList/" \
    -e "/FoamFile/,/}/s/funkyCellDecomp/cellDecomposition/" \
    -e "/dimensions/,/internalField/d" \
    -e "/boundaryField/,/^}/d" $fid > $fout

echo "Convert: $fid --> $fout"

total=`head -n50 $fout | grep "^(" -B1 | head -n1 `

echo "proc.(s): "
sed -n -e "/(/,/^;/p" $fout | sed -e '1d' -e '$d' | sed '$d' | \
sort -n | uniq -c | sed -e "s/^[ \t]*//" | while read line; do
    j=`echo $line | cut -f1 -d" "`
    count=`expr $count + $j`
    print_count $line $count `expr $total - $count`
done
return
}

function get_latest_time_dir()
{
   echo `ls -vd * | grep -E "^[0-9]" | sort -gr | head -n1`
return 
}

###############################################################3
# main body
###############################################################3

case "$1" in
    "-h"|"--help"|"-?"|"-help")
        usage;
        ;;
    *)
esac

if [ $# -ge 1 ]; then convert_to_labelList $1; exit; fi

create_dict

sed -n -e "/manualCoeffs/,/}/p" $fin | \
sed -e "s/^[ \t]*//" | sed -e "s/\/\/.*//" -e "s/;//" | \
grep "proc\." | while read line; do add_proc ${line#proc.}; done

execute_funky

DIR=$(get_latest_time_dir);
fid=`ls ${DIR}/funkyCellDecomp*`

convert_to_labelList $fid; echo


