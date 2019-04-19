#! /bin/bash

function usage()
{
cat << EOF

    Usage: `basename $0` [options] <sync.me>

    options:

        -h, --help      show this help message and exit
        -n, --dryrun    dryrun only i.e. no files are actually transferred
        -s, --show      show rsync command and exit
        -p, --print     create sample <sync.me> file and exit

EOF
}

function print_config()
{
cat > ${CFGFILE} << EOF
# sync.me (use hash symbol for comments)
# "rsync" options, please change as you see fit
# by defaults these options download from SERVER to current working directory

# please uncomment to modify
# OPTS="$OPTS"
# SERVER=$SERVER
# SSHKEY="$SSHKEY"
# SSHOPTS="$SSHOPTS"
USERNAME="$USERNAME"

# files to be downloaded defined in terms of SRC and DEST
# note:
#   - helper parameters can be created to simplify the definition of SRC and DEST
#   - by convention shell parameters are to be defined using CAPITAL LETTERS
#   - SRC are defined as a list (bash-syntax), hence the round bracket
#   - Use -s or -n option to show the rsync command
#       e.g.: download 2 files: SRC=(file1 file2) or using wildcard SRC=(file{1,2})
#       e.g.: more wildcard examples SRC=(folder*/log.run folder*/postProcessing)
ROOT="\${USERNAME}@\${SERVER}:/scratch/\${USERNAME}/add-your-src-folder-here"
SRC=(log.run*,postProcessing)
DEST="./"

EOF
cat $CFGFILE
}

CFGFILE="./sync.me"
OPTS="--recursive --times --progress --links --itemize-changes -s -R"
SERVER="liger.ec-nantes.fr"
USERNAME="<put-your-user-name-here>"
SSHKEY="~/.ssh/id_rsa"
SSHOPTS="-q"

POSITIONAL=()
if [ "$#" -eq 0 ]; then
    if [ ! -f "$CFGFILE" ]; then usage; exit; fi
fi
while [ "$#" -gt 0 ]; do
key=$1
case $key in
    -h|--help)
        usage; exit 0;
        shift
        ;;
    -n|--dryrun)
        OPTS+=" -n "
        SHOWCMD="yes"
        shift
        ;;
    -s|--show)
        SHOWCMD="yes"
        NOEXEC="yes"
        shift
        ;;
    -p|--print)
        print_config; exit 0
        shift
        ;;
    *)
        POSITIONAL+=("$key")
        shift
        ;;
esac
done
if [ "${#POSITIONAL[@]}" -eq 0 ]; then
    echo "Use default cfg.: ${CFGFILE}"
    POSITIONAL+=("$CFGFILE")
fi
set -- "${POSITIONAL[@]}"
if [ "$#" -ge 1 ]; then SYNCFILES=($@); fi

for iFile in "${SYNCFILES[@]}"; do
    if [ -f "${iFile}" ]; then
        source "${iFile}"
        ROOT=${ROOT%/}; ROOT=${ROOT%/.}; ROOT=${ROOT}/.
        for i in "${SRC[@]}"; do
            RUN="rsync ${OPTS} -e \"ssh $SSHOPTS -i $SSHKEY\" "
            RUN+="${ROOT}/${i} ${DEST}"
            if [ -n "$SHOWCMD" ]; then echo "bash: "$RUN; fi
            if [ -n "$NOEXEC" ]; then continue; fi
            eval $RUN
        done
    else
        echo "$(basename $0): file not found: ${iFile}"
    fi
done


