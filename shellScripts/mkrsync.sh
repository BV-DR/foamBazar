#! /bin/bash

if [ $# -lt 2 ]; then
cat << EOF

download file(s)/dir(s) from mek-cluster using rsync 

usage: $0 'remote_file(s)' destination

    remote              local
    -------------------------        
    new                 copy new file(s)
    deleted             ignored, don't delete local file(s)
    modified            copy/overwrite local file(s)

    symlink(s)          copy as symlink(s)
    timestamp           preserved
    directory           copy recursively

EOF
exit 1;
fi

HOSTIP=192.38.95.65
SRC=$1
DEST=$2

rsync -t -v -l -p -S -r sose@${HOSTIP}:${SRC} $DEST

