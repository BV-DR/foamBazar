#! /bin/bash

cat > .tmp << EOF
set title ""
set xlabel ""
set ylabel ""
plot "< cat run.log | grep Time | grep -v Ex | sed 's/Time = //'" using (\$1):(\$1) with lines title "time"
pause 1
reread
EOF

exec gnuplot .tmp -

