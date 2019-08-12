#!/bin/bash

#chmod u+x script.sh

# Execution
# $ bash execution.sh

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
rm master.zip
wget https://github.com/tcameronwaller/metabonet/archive/master.zip
unzip master.zip
cp -r ~/Downloads/metabonet-master/dock_template/ ~/
cd ~
mv dock_template dock_metabonet
ls "$path_dock/source/customization/"

##########
# Access source files.

echo "Recon2M.2"
cd ~/Downloads/
rm ./Recon2M.2_MNX_Entrez_Gene.xml
wget https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml
cp ~/Downloads/Recon2M.2_MNX_Entrez_Gene.xml $path_dock/source/recon2m2.xml
rm ./Recon2M.2_MNX_Entrez_Gene.xml

exit 0

echo "HMDB"
rm hmdb_metabolites.zip
rm hmdb_metabolites.xml
wget http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
unzip hmdb_metabolites.zip
cp ~/Downloads/hmdb_metabolites.xml $path_dock/source/hmdb_metabolites.xml
rm hmdb_metabolites.zip
rm hmdb_metabolites.xml

echo "measurements"

ls $path_dock/source/measurement/

##########
# Install MetaboNet


##########
# Curate metabolic model.

metabonet model -d $path_dock -r

echo "reconcile to MetaNetX"
ls $path_dock/reconciliation/

echo "curate and organize metabolic model"
ls $path_dock/source/customization/
metabonet model -d $path_dock -ceauv

echo "model files"
ls $path_dock/model/

##########
# Define metabolic networks.

echo "define metabolic networks"
ls $path_dock/source/customization/

echo "compartments true, hubs true"
metabonet network -d $path_dock -yc -np -v
path_network=$path_dock/network/compartments-true_hubs-true/
mkdir $path_network
mv $path_dock/network/links.pickle $path_network/
mv $path_dock/network/nodes_reactions.pickle $path_network/
mv $path_dock/network/nodes_metabolites.pickle $path_network/
mv $path_dock/network/network_cytoscape.json $path_network/
mv $path_dock/network/network_networkx.pickle $path_network/

echo "compartments true, hubs false"
metabonet network -d $path_dock -ycs -np -v
path_network=$path_dock/network/compartments-true_hubs-false/
mkdir $path_network
mv $path_dock/network/links.pickle $path_network/
mv $path_dock/network/nodes_reactions.pickle $path_network/
mv $path_dock/network/nodes_metabolites.pickle $path_network/
mv $path_dock/network/network_cytoscape.json $path_network/
mv $path_dock/network/network_networkx.pickle $path_network/



echo "compartments false, hubs true"
metabonet network -d $path_dock -y -np -v
path_network=$path_dock/network/compartments-false_hubs-true/
mkdir $path_network
mv $path_dock/network/links.pickle $path_network/
mv $path_dock/network/nodes_reactions.pickle $path_network/
mv $path_dock/network/nodes_metabolites.pickle $path_network/
mv $path_dock/network/network_cytoscape.json $path_network/
mv $path_dock/network/network_networkx.pickle $path_network/


echo "compartments false, hubs false"
metabonet network -d $path_dock -ys -np -v
path_network=$path_dock/network/compartments-false_hubs-false/
mkdir $path_network
mv $path_dock/network/links.pickle $path_network/
mv $path_dock/network/nodes_reactions.pickle $path_network/
mv $path_dock/network/nodes_metabolites.pickle $path_network/
mv $path_dock/network/network_cytoscape.json $path_network/
mv $path_dock/network/network_networkx.pickle $path_network/

echo "network files"
ls $path_dock/network/
