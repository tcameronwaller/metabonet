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

echo "model files"
ls $path_dock/model/

echo "measurement files"
ls $path_dock/measurement/

# Plot.

metabonet network -d $path_dock -t
