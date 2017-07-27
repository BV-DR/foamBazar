#! /bin/bash

function usage()
{
cat << EOF
Usage: `basename $0` src dest [time]
EOF
exit;
}

if [ $# -lt 2 ]; then usage; fi

SRC=$1; TARG=$2
if [ $# -eq 3 ]; then TIME=$3; fi
if [ -d "$TARG" ] ; then echo "target already existed: $TARG  ... abort"; usage; fi
if [ ! -d "$SRC" ]; then echo "src not found: $SRC  ... abort"; usage; fi

echo -e "Clone openfoam case\nfrom:" $SRC "\nto  :" $TARG
(
    mkdir -p $TARG
    cp -P "$SRC"/* "$TARG"
    cp -rP "$SRC"/{0,constant,system} "$TARG"
) > /dev/null 2>&1

if [ -d $SRC/processor0 ]; then
    echo "Clone processor*"
    for i in "$SRC"/processor*; do mkdir "$TARG"/`basename "$i"`; done
    (
        for i in "$SRC"/processor*; do
            cp -rpP "${i}"/{0,constant} "$TARG"/`basename "$i"`;
        done
    ) > /dev/null 2>&1
fi

if [ -n "$TIME" ]; then
    if [ -d "$SRC/$TIME" ]; then
        echo "Clone time (serial):" $TIME
        cp -rP "$SRC/$TIME" "$TARG"
    fi
    if [ -d "$SRC/processor0/$TIME" ]; then
        echo "Clone time (parallel): processor*/"$TIME
        for i in "$SRC"/processor*; do cp -rpP "${i}/$TIME" "$TARG"/`basename "$i"`; done
    fi
fi

rm -fr $TARG/log.*
touch $TARG/log.clone-`basename $SRC`

