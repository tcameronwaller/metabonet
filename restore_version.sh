#!/bin/bash

#chmod u+x script.sh

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous version of the program..."

cd ~/Downloads/
rm master.zip
rm -r metabonet-master/

# Access current version of the program.

echo "access current version of the program..."

wget https://github.com/tcameronwaller/metabonet/archive/master.zip
unzip master.zip
cd metabonet-master/
ls ./
sudo python3 setup.py install
cd ~
metabonet --help

