#! /bin/bash
if [ $# -lt 2 ]; then
cat << EOF

clone a foam case

Usage: $0 src-dir targ-dir

EOF
exit;
fi

SRC=$1; TARG=$2

if [ -d $TARG ]; then echo "target already existed: $TARG  ... abort"; exit; fi
if [ ! -d $SRC ]; then echo "src not found: $SRC  ... abort"; exit; fi

mkdir -p $TARG
cp -P $SRC/* $TARG 2>&1 > /dev/null
cp -rP $SRC/{0,constant,system} $TARG 2>&1 > /dev/null

if [ $# -eq 2 ]; then

if [ -d $SRC/processor0 ]; then
    for i in $SRC/processor*; do mkdir $TARG/`basename $i` 2>&1 > /dev/null; done
    for i in $SRC/processor*; do cp -rp ${i}/{0,constant} $TARG/`basename $i` 2>&1 > /dev/null; done
fi

fi

rm -fr $TARG/log.*
touch $TARG/log.clone-`basename $SRC`

