#! /bin/bash

if [ $# -lt 1 ]; then
cat << EOF

reset internalField

usage: $0 bc-file(s)

    volScalarField: reset with uniform 0;
    volVectorField: reset with uniform (0 0 0);

EOF
exit 1;
fi

for i in $*; do
    if [ ${i##*.} == "gz" ]; then
        rm -f ${i%.gz} > /dev/null 2>&1 
        gunzip $i > /dev/null 2>&1 
        i=${i%.gz}
    fi
    
    case `basename $i` in
    alpha1|gamma)
        if [ `grep internalField $i | grep "nonuniform List" | wc -l` == "1" ]; then
            sed -i -e "/internalField/,/;/d" $i
            sed -i -e "20 i internalField   uniform 0;" $i
        fi
        ;;
    *)
        echo "$0: vector or scalar? ... abort"
        exit 1; 
    esac
done

