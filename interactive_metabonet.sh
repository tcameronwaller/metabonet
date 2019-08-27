#!/bin/bash

#chmod u+x script.sh

# Interactive Execution Script
# $ bash execution.sh

#######################
# Get current metabonet
# directory from user
#######################
echo "+------------------------------------------------------------------------------+" 2>/dev/null
echo "|                                  MetaboNet                                   |"
echo "+------------------------------------------------------------------------------+" 2>/dev/null
echo "| Copyright (C) 2019 Thomas Cameron Waller                                     |"
echo "| This project supports curation of the model of human cellular metabolism     |"
echo "| definition of custom metabolic networks to represent this model, and         |"
echo "| analyses of these networks.                                                  |"
echo "|                                                                              |" 2>/dev/null
echo "| If you use Metabonet, please cite the following:                             |"
echo "| T. Cameron Waller, Jordan A. Berg, Brian E. Chapman, Jared P. Rutter. 2019.  |"
echo "| \"Compartments and Hubs Differentiate the Relevance of Metabolic Networks to  |"
echo "| Metabolomic Experiments\".                                                    |"
echo "|                                                                              |"
echo "+------------------------------------------------------------------------------+" 2>/dev/null
echo "| This interactive script will walk you through installing MetaboNet curating  |"
echo "| modeling the human metabolic network, and analyzing data using this model.   |"
echo "|                                                                              |"
echo "| Note: Requires PyPi (comes with most current distributions of Python)        |"
echo "|                                                                              |"
echo "+------------------------------------------------------------------------------+" 2>/dev/null
echo "|                                                                              |" 2>/dev/null
echo "+------------------------------------------------------------------------------+" 2>/dev/null
echo "| Step 1: Download Source Files"
echo "| Here is your current working directory:" 2>/dev/null
pwd 2>/dev/null
echo "| " 2>/dev/null
read -p "| Are you in the MetaboNet directory? (Y/n)" -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
  METABONET="."
else
  echo ""
  read -p "| Please provide the full path to the MetaboNet directory: " -r
  METABONET=$REPLY
fi

# Check for path formatting
if [[ ${METABONET: -1} == "/" ]]
then
  METABONET=${METABONET::-1}
else
  METABONET=$METABONET
fi
echo "" 2>/dev/null
echo "+-------------------------------------------------------------------------------" 2>/dev/null

# Set dock path
DOCK="$METABONET/dock_template"

#######################
# Access source files
#######################
RECON='https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml'

rm -rf $DOCK/source/recon2m2.xml 2>/dev/null
rm -rf $DOCK/source/Recon2M.2_MNX_Entrez_Gene.xml 2>/dev/null

echo "| Currently set to download the Recon database from: "
echo "| $RECON"
read -p "| Do you want to download and use this version? (Y/n)" -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
  RECON=$RECON
else
  echo "|"
  read -p "| Please provide the URL for the version you would like to use: " -r
  RECON=$REPLY
fi
echo ""
echo "| Downloading Recon..."
curl -L $RECON -o $DOCK/source/recon2m2.xml


HMDB='http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'

rm -rf $DOCK/source/hmdb_metabolites.zip
rm -rf $DOCK/source/hmdb_metabolites.xml

echo "| Currently set to download and use the HMDB database from: "
echo "| $HMDB"
read -p "| Do you want to download this version? (Y/n)" -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
  HMDB=$HMDB
else
  echo "|"
  read -p "| Please provide the URL for the version you would like to use: " -r
  HMDB=$REPLY
fi
echo ""
echo "| Downloading HMDB..."
curl -L $HMDB -o $DOCK/source/hmdb_metabolites.zip
unzip $DOCK/source/hmdb_metabolites.zip -d $DOCK/source/
rm -rf $DOCK/source/hmdb_metabolites.zip

echo "| "
echo "| Source files downloaded."
echo "|"
echo "+-------------------------------------------------------------------------------" 2>/dev/null


#######################
# Check user has
# included all desired
# measurement records
#######################
echo "| Step 2: Include Measurement Records"
echo "| Current measurement records include:"
ls $DOCK/source/measurement/
echo "| Check that no desired measurement records are missing, then press any key."
read -p  "| Include any missing files or folders at $DOCK/source/measurements" -n 1 -r
echo ""
echo "|"
echo "+-------------------------------------------------------------------------------" 2>/dev/null

#######################
# Install MetaboNet
#######################
echo "| Step 3: Install MetaboNet"
echo "| Installing MetaboNet..."
python $METABONET/setup.py install
echo "| MetaboNet installed."
echo "|"
echo "+-------------------------------------------------------------------------------" 2>/dev/null


#######################
# Curate metabolic
# model
#######################
echo "| Step 4: Curate Metabolic Model"
echo "| Curating MetaboNet metabolic model..."

mkdir $DOCK/reconciliation/
METANETX="$DOCK/source/reconciliation_2019-08-13"
echo "| Currently reconciling network using archived MetaNetX network file: "
echo "| $METANETX"
read -p "| Do you wish to run a new reconciliation? (Y/n)" -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]
then
  echo ""
  echo "| Preparing reconciliation..."
  metabonet model -d $DOCK -r
  echo ""
  echo "| The following will walk you through custom reconciliation of the metabolic "
  echo "| network. We will provide you with each task step-by-step. Press any key to "
  echo "| advance to the next step."
  echo "| "
  echo "| - Import the file $DOCK/reconciliation/recon2m2_reconciliation.xml "
  read -p "| into MetaNetX @ https://www.metanetx.org/ and run. (press any key to continue)" -n 1 -r
  echo ""
  echo "| "
  echo "| - Save _Mapping summary_ from [MetaNetX][15] as the following file."
  echo "| Save the mapping summary to $DOCK/reconciliation/metanetx_import_report.tsv "
  echo "| This summary is useful for review and also includes information to curate "
  read -p "| (press any key to continue)" -n 1 -r
  echo ""
  echo "| "
  echo "| - Save files from MetaNetX to the following directories and files:"
  echo "| (press any key to continue for each file)"
  read -p "| * $DOCK/reconciliation/recon2m2_metanetx_compartments.tsv " -n 1 -r
  echo ""
  read -p "| * $DOCK/reconciliation/recon2m2_metanetx_genes.tsv " -n 1 -r
  echo ""
  read -p "| * $DOCK/reconciliation/recon2m2_metanetx_metabolites.tsv " -n 1 -r
  echo ""
  read -p "| * $DOCK/reconciliation/recon2m2_metanetx_reactions.tsv " -n 1 -r
  echo ""
else
  cp -r $DOCK/source/reconciliation_2019-08-13/ $DOCK/
  mv $DOCK/reconciliation_2019-08-13/ $DOCK/reconciliation/
  ls $DOCK/reconciliation/
fi
METANETX="$DOCK/reconciliation"
echo "| "
echo "| Customization files:"
ls $DOCK/source/customization/
echo "| "
echo "| If you would like to edit any customizable files, do so now @ $DOCK/source/customization/"
read -p "| (press any key to continue)" -n 1 -r
echo ""
echo "| "
echo "| Curating and organizing metabolic model..."

metabonet model -d $DOCK -ceauv

echo "| "
echo "| Model files:"
ls $DOCK/model/

echo "| "
echo "| Curating and organizing measurements..."
metabonet model -d $DOCK -m

echo "| "
echo "| Measurement files:"
ls $DOCK/measurement/

#######################
# Define metabolic
# networks
#######################
echo "| "
echo "| Defined metabolic networks:"
ls $DOCK/source/customization/

#######################
# Build metabolic
# model
#######################
echo "+-------------------------------------------------------------------------------" 2>/dev/null
echo "| Step 5: Build Metabolic Model "
echo "| We will now build the metabolic model."
echo "| "
echo "| Building model..."

echo "| compartments true, hubs true "
metabonet network -d $DOCK -yc -np -v
echo "| "

echo "| Compartments true, Hubs false "
metabonet network -d $DOCK -ycs -np -v
echo "| "

echo "| Compartments false, Hubs true "
metabonet network -d $DOCK -y -np -v
echo "| "

echo "| Compartments false, Hubs false "
metabonet network -d $DOCK -ys -np -v
echo "| "
echo "+-------------------------------------------------------------------------------" 2>/dev/null
echo "| Step 6: Integrate measurements into model"
echo "| Integrating measurements to network's metabolites..."
metabonet network -d $DOCK -m


#######################
# Analysis
#######################
echo "| "
echo "+-------------------------------------------------------------------------------" 2>/dev/null
echo "| Step 7: Analyzing Network"
echo "| Analyzing network..."
metabonet network -d $DOCK -a
mv $DOCK/network $path_network
#mv $DOCK/network/links.pickle $path_network/
#mv $DOCK/network/nodes_reactions.pickle $path_network/
#mv $DOCK/network/nodes_metabolites.pickle $path_network/
#mv $DOCK/network/network_cytoscape.json $path_network/
#mv $DOCK/network/network_networkx.pickle $path_network/
#mv $DOCK/network/measurement/ $path_network/
#mv $DOCK/network/analysis/ $path_network/
echo "| "
echo "| Network files:"
ls $DOCK/

# Plot
echo "| "
echo "| Plotting results..."
metabonet network -d $DOCK -t
echo "| "
echo "+------------------------------------------------------------------------------+" 2>/dev/null
echo "|                         MetaboNet execution complete                         |"
echo "+------------------------------------------------------------------------------+" 2>/dev/null
