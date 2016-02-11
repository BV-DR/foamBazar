#! /bin/bash

# Name of the script
MY_NAME=`basename $0`;

function usage ()
{
    if [ -n "$1" ]; then
        echo "Unknown option : $1"
    fi

cat << EOF

Usage: $MY_NAME [option(s)] imagefile(s)

    options:
    
    -dt=deltaT      dt pr. frame, default: 0.04 (25 fps)
    -tdir=/path/    path to time dir
    -i              don't create dir.timestamp, overwrite orig. (edit in place)

EOF
}

function format_hh_mm_ss_sss ()
{
    hh=`echo $1 | awk '{ printf("%02.f",$1/3600) }'`
    mm=`echo $1 $hh | awk '{ printf("%02.f",$1/60-60*$2) }'`
    ss=`echo $1 $hh $mm | awk '{ printf("%02.f",$1-3600*$2-60*$3) }'`
    sss=${1#*.};
    echo "$hh:$mm:$ss.$sss"
}

deltaT="0.04"; tdir=
overwrite=0; files=
suffix='.timestamp'
opts='-pointsize 32 -fill white -verbose'

while [ -n "$1" ]; do
    case $1 in
        -dt=*)
            deltaT=`echo $1 | cut -d '=' -f2-`
            shift
            ;;
        -tdir=*)
            tDir=`echo $1 | cut -d '=' -f2-`
            shift
            ;;
        -i)
            overwrite=1; prefix="";
            shift
            ;;
        -*)
            usage $1
            exit 1
            ;;
        *)
            files+="$1 "
            shift
            ;;
    esac
done


for i in $files; do

    frame=${i#*.}; frame=`echo ${frame%.*} | bc`; text=`echo "scale=${#deltaT}; $frame * $deltaT" | bc -l`
    hh=`echo $text | awk '{ printf("%02f",$1/3600) }'`; hh=`echo ${hh%.*} | awk '{printf("%02d",$1)}'`
    mm=`echo $text $hh | awk '{ printf("%02f",$1/60-60*$2) }'`; mm=`echo ${mm%.*} | awk '{printf("%02d",$1)}'`;
    ss=`echo $text $hh $mm | awk '{ printf("%02f",$1-3600*$2-60*$3) }'`; ss=`echo ${ss%.*} | awk '{printf("%02d",$1)}'`
    sss=${text#*.};
    convert $i $opts -annotate +8+32 "$hh:$mm:$ss.$sss, $frame" $i${suffix}.${i##*.}
done



#echo $files

exit

if [ $# -lt 1 ]; then
cat << EOF
Usage: $0 jpg.prefix fps vcodec

    jpg.prefix      for now it is a directory, will be fixed soon
    fps             default 25 fps
    vcodec          default huffyuv, see mencoder man page(s)

EOF
exit;
fi

if [ $# -lt 1 ]; then
cat << EOF
Usage: $0 jpg.prefix fps vcodec

    jpg.prefix      for now it is a directory, will be fixed soon
    fps             default 25 fps
    vcodec          default huffyuv, see mencoder man page(s)

EOF
exit;
fi



exit
for file in *.jpg; do
convert $file -font arial.ttf -fill -annotate +100+100 new-$file
done

convert pd/pd.0001.jpg -pointsize 32 -fill white -annotate "+32+64" "00:00:00.003" 1.jpg


