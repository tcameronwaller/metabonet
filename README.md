# [MetaboNet](1)

This project supports curation of the human metabolic models, definition of
custom metabolic networks, and analyses of these networks.

## Launch MetaboNet's [Singularity](2) Container

```bash
docker pull jordanberg/metabonet:publication
docker run metabonet
```

## Install MetaboNet

Download [MetaboNet](1) from [GitHub](3).

```bash
$ cd ~/Downloads/
$ wget https://github.com/tcameronwaller/metabonet/archive/master.zip
$ unzip master.zip
$ cd metabonet-master/
```

Install MetaboNet to default directory for third-party Python packages.

```bash
$ sudo python3 setup.py install
or
$ pip3 install -e setup.py
```

Confirm installation. If no errors are output, installation was successful.

```bash
$ metabonet --help
$ metabonet model --help
$ metabonet network --help
$ metabonet clean --help
```

## Set Up Dock Directory

Procedures in MetaboNet utilize a central directory in which to organize both
source and product files. These instructions refer to this directory by the
name "dock". This directory can be anywhere the user wants.

A template "dock" directory with default files is accessible within the
MetaboNet repository.

## Set Up Source Information Within Dock

### Metabolic Information

MetaboNet derives information from both a model of human metabolism and from a
database of information about metabolites. As the necessary files are large,
these instructions specify how to access them from original repositories.

Download model of human metabolism.

MetaboNet derives information about human metabolism from model
[Recon 2M.2](4). Access [Recon 2M.2](4) by downloading file
"Recon2M.2_MNX_Entrez_Gene.xml" (14.2 MB) from record [583326](5) on the
[Zenodo](6) repository. Change the name of this file to "recon2m2.xml" and
place it in the following location.

```bash
$ wget https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml?download=1
$ mv Recon2M.2_MNX_Entrez_Gene.xml?download=1 recon2m2.xml
$ mv recon2m2.xml dock/source/
```

Download database about human metabolites.

MetaboNet enhances information about human metabolites from the
[Human Metabolome Database](7) [HMDB](8). Access the [HMDB](8) by
[downloading](9) version "4.0" of file "hmdb_metabolites.xml" (4.2 GB).
Decompress this file and place it in the following location.

```bash
$ wget http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
$ unzip hmdb_metabolites.zip
$ mv hmdb_metabolites.xml dock/source/
$ rm hmdb_metabolites.zip
```

### Customization Parameters

MetaboNet is versatile in its support for customization at many levels. Files
of default customization parameters are available within the "customization"
directory in the [MetaboNet](1) repository on [GitHub](3).

Copy customization files to "dock" directory.

```bash
$ cp ~/Downloads/metabonet-master/customization/ /dock/source/
```

### Metabolomic Measurements

MetaboNet curates and analyzes metabolomic measurements. The original data sets
come from [Metabolomics Workbench](10), a [repository](11) for data from metabolomic
studies. The user has 2 options to prepare these data sets for use in MetaboNet.

#### 1. Access and prepare metabolomic measurements

1. Locate specific project and study on [Metabolomics Workbench][11].
2. Copy information about samples.
-Click on study's tab "Show all samples".
-Copy information from sample table.
-Save file as "samples.tsv".
3. Copy information about analytes.
-Click on study's tab "Show named metabolites".
-Save file as "analytes.tsv".
4. Copy information about measurements.
-Click on study's tab "Download named metabolite data".
-Save file as "measurements.tsv".
5. Copy information about signals.
-Click on study's tab "Download all metabolite data".
-Save file as "signals.tsv".
6. Correct errors.
Project PR000305, Study ST000390
analytes.tsv
"ascorbic acid", "Vitamin C", change PubChem identifier from 5785 to 54670067.
Project PR000322, Study ST000412
analytes.tsv
"putrescine", "Putrescine", change PubChem identifier from 1049 to 1045.
Remove record for "maltose", "Maltose", PubChem 439186.
That record is redundant with the record for "cellobiose", PubChem 6255.

#### 2. Copy measurement files to "dock" directory

```bash
$ cp ~/Downloads/metabonet-master/measurement/ /dock/source/
```

## Curate Metabolic Model and Metabolomic Measurements

MetaboNet's "model" routine includes multiple procedures to curate the
metabolic model and to curate metabolomic measurements. This curation prepares
to define custom metabolic networks and integrate measurements with these
networks.

### Reconcile Metabolic Model for Import to MetaNetX

Execute the "reconciliation" procedure of the "model" routine in MetaboNet.

```bash
$ metabonet model -d /dock/ -r
```

[MetaNetX](12) is a [repository](13) of metabolic models with tools to curate
these models.

1. Import the following file to [MetaNetX](13).
```bash
/dock/reconciliation/recon2m2_reconciliation.xml
```
2. Save "Mapping summary" as the following file. This summary is useful for
review and also includes useful information to curate names of metabolites.
Import category "chemical compounds mapped to the MNXref namespace with
different ID and different description" is of particular interest for curation
of information about metabolites.
```bash
/dock/reconciliation/metanetx_import_report.tsv
```
3. Save files from MetaNetX to the following directories and files.
```bash
compartments.tsv -> /dock/reconciliation/recon2m2_metanetx_compartments.tsv
enzymes.tsv -> /dock/reconciliation/recon2m2_metanetx_genes.tsv
chemicals.tsv -> /dock/reconciliation/recon2m2_metanetx_metabolites.tsv
reactions.tsv -> /dock/reconciliation/recon2m2_metanetx_reactions.tsv
```

### Curate Metabolic Model

Execute the "collection", "extraction", "enhancement", "curation", and
"conversion" procedures of the "model" routine in MetaboNet.

The "enhancement" procedure alone requires about 1 hour to complete.

```bash
$ metabonet model -d /dock/ -ceauv
```

A file of special use for automatic curation of reactions follows.

```bash
/dock/enhancement/reactions_filter.tsv
```

Customize curation by editing the following files of parameters.

```bash
/dock/source/customization/curation_compartments.tsv
/dock/source/customization/curation_processes.tsv
/dock/source/customization/curation_reactions.tsv
/dock/source/customization/curation_metabolites.tsv
```

### Curate Metabolomic Measurements

Execute the "measurement" procedure of the "model" routine in MetaboNet.

```bash
$ metabonet model -d /dock/ -m
```

Review the following files to check the matches between analytes and
metabolites. Discrepancies might justify modifications to curation of
metabolites.

```bash
/root/measurement/study_one_report.tsv
/root/measurement/study_two_report.tsv
/root/measurement/study_three_report.tsv
/root/measurement/study_four_report.tsv
/root/measurement/study_five_report.tsv
```

## Define Custom Metabolic Networks

Customize definition of metabolic networks by editing the following files of
parameters.

```bash
/dock/source/customization/filtration_compartments.tsv
/dock/source/customization/filtration_processes.tsv
/dock/source/customization/simplification_reactions.tsv
/dock/source/customization/simplification_metabolites.tsv
```

Execute the "candidacy", "network", and "conversion" procedures of the
"network" routine in MetaboNet.

compartments true, hubs true
```bash
$ metabonet network -d /dock/ -yc -np -v
```

compartments true, hubs false
```bash
$ metabonet network -d /dock/ -ycs -np -v
```

compartments false, hubs true
```bash
$ metabonet network -d /dock/ -y -np -v
```

compartments false, hubs false
```bash
$ metabonet network -d /dock/ -ys -np -v
```

## Analyze Metabolic Networks

Execute the "analysis" procedure of the "network" routine in MetaboNet.

```bash
$ metabonet network -d /dock/ -a
```

Run time of "analysis" procedure can be 0.5-2.0 hours depending on the network.

## Visualize Metabolic Networks

```bash
/dock/conversion/network_elements_cytoscape.json
```

## Integrate Metabolomic Measurements in Metabolic Network

Execute the "measurement" proceduree of the "network" routine in MetaboNet.

```bash
$ metabonet network -d /dock/ -m
```

The following file includes measurements from all studies ready for integration
with the metabolic network in Cytoscape.

```bash
/root/measurement/metabolites.tsv
```

## Removal of automatically-generate files and directories

```bash
$ metabonet clean -d /dock/
```





[1]: [https://github.com/tcameronwaller/metabonet]
[2]: [https://www.sylabs.io/docs/]
[3]: [https://github.com/]
[4]: [https://www.ncbi.nlm.nih.gov/pubmed/29078384]
[5]: [https://zenodo.org/record/583326]
[6]: [https://zenodo.org/]
[7]: [https://www.ncbi.nlm.nih.gov/pubmed/29140435]
[8]: [http://www.hmdb.ca/]
[9]: [http://www.hmdb.ca/downloads]
[10]: [https://www.ncbi.nlm.nih.gov/pubmed/26467476]
[11]: [http://www.metabolomicsworkbench.org/]
[12]: [https://www.ncbi.nlm.nih.gov/pubmed/26527720]
[13]: [https://www.metanetx.org/]
