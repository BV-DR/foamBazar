#! /bin/bash
if [ $# -lt 1 ]; then
cat << EOF
Usage: $0 jpg.prefix fps vcodec

    jpg.prefix      for now it is a directory, will be fixed soon
    fps             default 25 fps
    vcodec          default huffyuv, see mencoder man page(s)

EOF
exit;
fi


if [ $# -lt 2 ]; then FPS=25; else FPS=$2; fi
if [ $# -lt 3 ]; then VCODEC=huffyuv; else VCODEC=$3; fi
echo "FPS: $FPS"
echo "VCODEC: $VCODEC"
mencoder "mf://${1}*.png" -mf type=png:fps=$FPS -o ${1}.avi -nosound -ovc lavc -lavcopts vcodec=$VCODEC:format=422p

# mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -nosound -o output.avi

# mencoder "mf://./${1}*.png" -mf type=png:fps=10 -o ./output.avi -nosound -ovc lavc -lavcopts vcodec=huffyuv:format=422p

