#!/usr/bin/env python

import os, sys
import argparse
import configparser
import subprocess as sub

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', dest='inputFile', default='sync.me', help='read source directories from input file.')
parser.add_argument('-u','--update', dest='updateTree', action='store_true', help='use this option to update list of directories in tree.')
parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print sample config file')
args = parser.parse_args()

#Print config file
if args.showConfig:
    print('Output default parameters to file: "sync.me"')
    with open('sync.me','w') as f:
        f.write('[sync]\n')
        f.write('server = username@liger.ec-nantes.fr\n')
        f.write('sshkey = ~/.ssh/id_rsa\n') 
        f.write('src = /scratch/username/mycase\n')
    os._exit(1)

#Read inputfile
config = configparser.ConfigParser()
config.read(args.inputFile)
server = str(config['sync']['server'])
sshkey = str(config['sync']['sshkey'])
src = str(config['sync']['src'])

#Update tree of folders to download
if args.updateTree or (not os.path.isfile('.filetree')):
    print('Read and update folders to synchronize')
    with open('.filetree','w') as ftree:
        sub.call('ssh -a '+server+' "cd '+src+'; find . -type d -name postProcessing"', stdout=ftree, shell=True)

#Read folders to download
with open('.filetree') as f:
    tree = f.readlines()
tree = [x.strip()[2:] for x in tree if x.startswith('./')]

#Make all directories
print('Synchronize directories')
for f in tree:
    directory = os.path.dirname(f)
    if len(directory)>0:
        if not (os.path.exists(directory)): os.makedirs(directory)
        string = '( cd '+directory+'; scp -r '+server+':'+os.path.join(src,f)+' .; )'
        string = '( cd '+directory+'; rsync --recursive --times --progress --links --backup --itemize-changes -s -e "ssh -a -q -i '+sshkey+'" '+server+':'+os.path.join(src,f)+' ./ ;)'
    else:
        string = '(rsync --recursive --times --progress --links --backup --itemize-changes -s -e "ssh -a -q -i '+sshkey+'" '+server+':'+os.path.join(src,f)+' ./ ;)'
    sub.call(string,shell=True)
