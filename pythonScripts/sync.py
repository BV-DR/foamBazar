#!/usr/bin/env python

import os, sys
import argparse
import configparser
import subprocess as sub

server = "pepppewdd@liger.ec-nantes.fr"

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='inputFile', default='sync.me', help='read source directories from input file.')
parser.add_argument('-u','--update', dest='updateTree', action='store_true', help='use this option to update list of directories in tree.')
parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print sample config file')
args = parser.parse_args()

#Print config file
if args.showConfig:
    print('\nOutput default parameters to file: "sync.me"')
    with open('sync.me','w') as f:
        f.write('[sync]\n')
        f.write('server = username@liger.ec-nantes.fr\n')
        f.write('src = /scratch/username/mycase\n')
    os._exit(1)

#Read inputfile
config = configparser.ConfigParser()
config.read(args.inputFile)
server = str(config['sync']['server'])
src = str(config['sync']['src'])

#Update tree of folders to download
if args.updateTree or (not os.path.isfile('.filetree')):
    with open('.filetree','w') as ftree:
        sub.call('ssh '+server+' "cd '+src+'; find . -type d -name postProcessing"', stdout=ftree, shell=True)

#Read folders to download
with open('.filetree') as f:
    tree = f.readlines()
tree = [x.strip()[2:] for x in tree if x.startswith('./')]

#Make all directories
for f in tree:
    directory = os.path.dirname(f)
    if not os.path.exists(directory): os.makedirs(directory)
    string = '( cd '+directory+'; scp -r '+server+':'+os.path.join(src,f)+' .; )'
    #string = '( cd '+directory+'; rsync --recursive --times --progress --links --backup --itemize-changes -s -e ssh -q -i /root/.ssh/id_rsa '+server+':'+os.path.join(src,f)+' ./ ;)'
    sub.call(string,shell=True)



#CMD="--recursive --times --progress --links --backup --itemize-changes -s "
#for iFile in "${SYNCFILES[@]}"; do
#    if [ -f "${iFile}" ]; then
#        source "${iFile}"
#        for i in "${SRC[@]}"; do
#            rsync ${CMD} -e "ssh -q -i /root/.ssh/id_rsa" ${i} ${DEST}
#        done
#    else
#        echo "$(basename $0): file not found: ${iFile}"
#    fi
#done
