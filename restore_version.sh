#!/bin/bash

#chmod u+x script.sh

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous version of the program..."

cd ~/Downloads/
rm master.zip
#chmod -R 777 metabonet-master/
#rm -rf metabonet-master/

# Access current version of the program.

echo "access current version of the program..."

rm master.zip
wget https://github.com/tcameronwaller/metabonet/archive/master.zip
unzip master.zip
sudo chmod -R 777 ~/Downloads/metabonet-master/
cd metabonet-master/
ls ./

echo "................................................."
echo "................................................."
echo "................................................."
echo "................................................."
echo "................................................."
echo "uninstall previous version..."
python3 setup.py install --user --record installation_files.txt
cat installation_files.txt | sudo xargs rm -rf
rm -rf installation_files.txt

python3 setup.py install --user --force --prefix=
cd ~/Downloads/
rm master.zip
chmod -R 777 ~/Downloads/metabonet-master/
rm -rf metabonet-master/

cd ~
metabonet --help




