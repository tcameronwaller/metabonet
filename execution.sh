#!/bin/bash

#chmod u+x script.sh

# Execution
# $ bash execution.sh
# nohup bash execution.sh > ~/metabonet_report.txt &

path_dock=~/dock_metabonet

# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "Here are your working directories."
echo "dock: $path_dock"
echo "----------"
echo "--------------------------------------------------"

##########
# Echo each command to console.
set -x

echo "Organize source files and default parameters."

rm -rf $path_dock

#mkdir "$path_dock"

##########
# Access default parameters from repository

echo "parameters"
cd ~/Downloads/
rm -rf ~/Downloads/master.zip
rm -rf ~/Downloads/metabonet-master
wget https://github.com/tcameronwaller/metabonet/archive/master.zip
unzip master.zip
cp -r ~/Downloads/metabonet-master/dock_template/ ~/
cd ~
mv dock_template dock_metabonet
ls "$path_dock/source/customization/"

rm -rf ~/Downloads/master.zip
rm -rf ~/Downloads/metabonet-master

##########
# Access source files.

echo "Recon2M.2"
cd ~/Downloads/
rm ./Recon2M.2_MNX_Entrez_Gene.xml
wget https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml
cp ~/Downloads/Recon2M.2_MNX_Entrez_Gene.xml $path_dock/source/recon2m2.xml
rm ~/Downloads/Recon2M.2_MNX_Entrez_Gene.xml

echo "HMDB"
rm hmdb_metabolites.zip
rm hmdb_metabolites.xml
wget http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
unzip ~/Downloads/hmdb_metabolites.zip
cp ~/Downloads/hmdb_metabolites.xml $path_dock/source/hmdb_metabolites.xml
rm ~/Downloads/hmdb_metabolites.zip
rm ~/Downloads/hmdb_metabolites.xml

echo "measurements"

ls $path_dock/source/measurement/

##########
# Install MetaboNet



##########
# Curate metabolic model.

#metabonet model -d $path_dock -r

echo "reconcile to MetaNetX"
cp -r $path_dock/source/reconciliation_2019-08-13/ $path_dock/
mv $path_dock/reconciliation_2019-08-13/ $path_dock/reconciliation/
ls $path_dock/reconciliation/

echo "curate and organize metabolic model"
ls $path_dock/source/customization/
metabonet model -d $path_dock -ceauv

echo "model files"
ls $path_dock/model/

echo "curate and organize measurements"
metabonet model -d $path_dock -m

echo "measurement files"
ls $path_dock/measurement/

##########
# Define metabolic networks.

echo "define metabolic networks"
ls $path_dock/source/customization/

echo "compartments true, hubs true"
metabonet network -d $path_dock -yc -np -v
echo "integrate measurements to network's metabolites"
metabonet network -d $path_dock -m
echo "analyze network"
metabonet network -d $path_dock -a
path_network=$path_dock/compartments-true_hubs-true
mv $path_dock/network $path_network
#mv $path_dock/network/links.pickle $path_network/
#mv $path_dock/network/nodes_reactions.pickle $path_network/
#mv $path_dock/network/nodes_metabolites.pickle $path_network/
#mv $path_dock/network/network_cytoscape.json $path_network/
#mv $path_dock/network/network_networkx.pickle $path_network/
#mv $path_dock/network/measurement/ $path_network/
#mv $path_dock/network/analysis/ $path_network/

echo "compartments true, hubs false"
metabonet network -d $path_dock -ycs -np -v
echo "integrate measurements to network's metabolites"
metabonet network -d $path_dock -m
echo "analyze network"
metabonet network -d $path_dock -a
path_network=$path_dock/compartments-true_hubs-false
mv $path_dock/network $path_network

echo "compartments false, hubs true"
metabonet network -d $path_dock -y -np -v
echo "integrate measurements to network's metabolites"
metabonet network -d $path_dock -m
echo "analyze network"
metabonet network -d $path_dock -a
path_network=$path_dock/compartments-false_hubs-true
mv $path_dock/network $path_network

echo "compartments false, hubs false"
metabonet network -d $path_dock -ys -np -v
echo "integrate measurements to network's metabolites"
metabonet network -d $path_dock -m
echo "analyze network"
metabonet network -d $path_dock -a
path_network=$path_dock/compartments-false_hubs-false
mv $path_dock/network $path_network

echo "network files"
ls $path_dock/

# Plot.

metabonet network -d $path_dock -t
