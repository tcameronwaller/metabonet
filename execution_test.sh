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

#mkdir "$path_dock"

##########
# Access default parameters from repository


##########
# Access source files.


##########
# Install MetaboNet



##########
# Curate metabolic model.

##########
# Define metabolic networks.

echo "define metabolic networks"
ls $path_dock/source/customization/

echo "compartments false, hubs false"
metabonet network -d $path_dock -ys -np -v
echo "integrate measurements to network's metabolites"
metabonet network -d $path_dock -m
echo "analyze network"
metabonet network -d $path_dock -a
path_network=$path_dock/compartments-false_hubs-false
mv $path_dock/network $path_network
