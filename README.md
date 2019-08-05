# [MetaboNet][1]

This project supports curation of the model of human cellular metabolism,
definition of custom metabolic networks to represent this model, and analyses
of these networks.



## Install MetaboNet

**_DO:_** Download [MetaboNet][1] from [GitHub][2] and extract the archive.

```bash
$ cd ~/Downloads/
$ wget https://github.com/tcameronwaller/metabonet/archive/master.zip
$ unzip master.zip
$ cd metabonet-master/
$ ls ./
dock LICENSE metabonet README.md setup.py
```

**_DO:_** Inspect the repository's contents.

```bash
$ ls ~/Downloads/metabonet-master/
dock LICENSE metabonet README.md setup.py
```

**_DO:_** Install MetaboNet to default directory for third-party Python
packages.

```bash
$ cd metabonet-master/
$ sudo python setup.py install
or
$ pip install -e setup.py
```

**_DO:_** Confirm installation.

```bash
$ metabonet --help
$ metabonet model --help
$ metabonet network --help
$ metabonet clean --help
```



## Set Up Dock Directory

Procedures in [MetaboNet][1] organize both source (input) and product (output)
files within a central directory. These instructions refer to this directory by
the name _"dock"_. The user can place this directory wherever on their machine
but must tell [MetaboNet][1] where to find it.

A template _"dock"_ directory with default files is accessible within the
[MetaboNet][1] repository.

**_DO:_** Copy default _"dock"_ directory to accessible location.

```bash
$ cp ~/Downloads/metabonet-master/dock/ ~/dock/
$ ls ~/dock/
```

### Customizable Parameters

MetaboNet supports customization at the level of curation of the metabolic
model and at the level of definition of metabolic networks.

The [MetaboNet][1] repository on [GitHub][2] contains default parameters within
the "customization" directory.

**_DO:_** Inspect these files of customizable parameters.

```bash
$ ls ~/dock/source/customization/
```

### Metabolic Information

[MetaboNet][1] derives information from the [Recon 2M.2][3] model of human
metabolism and from the [Human Metabolome Database][6] [(HMDB)][7] of
metabolites. As the necessary files are large, these instructions specify how
to access them from the original repositories.

Access [Recon 2M.2][3] in file _"Recon2M.2_MNX_Entrez_Gene.xml"_
(14.2 MB) from record [583326][4] on the [Zenodo][5] repository.

**_DO:_** Download [Recon 2M.2][3] model of human metabolism.

```bash
$ cd ~/Downloads/
$ wget https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml
```

**_DO:_** Change the name of this file to _"recon2m2.xml"_ and move it to the
_"dock"_ directory.

```bash
$ mv ~/Downloads/Recon2M.2_MNX_Entrez_Gene.xml ~/dock/source/recon2m2.xml
$ ls ~/dock/source/
```

Access the [HMDB][6] in version _"4.0"_ of file _"hmdb_metabolites.xml"_
(4.2 GB).

**_DO:_** Download the [HMDB][6] database about human metabolites. Extract this
file and move it to the _"dock"_ directory.

```bash
$ cd ~/Downloads/
$ wget http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
$ unzip hmdb_metabolites.zip
$ mv ~/Downloads/hmdb_metabolites.xml ~/dock/source/hmdb_metabolites.xml
$ ls ~/dock/source/
```

### Metabolomic Measurements

[MetaboNet][1] integrates metabolomic measurements with metabolic networks. The
authors curated metabolomic measurements from five specific studies and used
these as examples. The original data sets for these metabolomic measurements
came from [Metabolomics Workbench][8], a [repository][9] for data from
metabolomic studies. These instructions explain the collection and curation of
these data sets for use in [MetaboNet][1], and they also explain how to
access ready archives.

#### Option 1: Access measurement files in _"dock"_ directory

```bash
$ ls ~/Downloads/dock/source/measurement/
```

#### Option 2: Access and prepare metabolomic measurements

1. Locate specific projects and studies on [Metabolomics Workbench][9].
- Project: [PR000305][10], Study: ST000390
- Project: [PR000058][11], Study: ST000061
- Project: [PR000322][12], Study: ST000412
- Project: [PR000322][12], Study: ST000412
- Project: [PR000599][13], Study: ST000842
2. Copy information about samples.
- Click on study's tab "Show all samples".
- Copy information from sample table.
- Save file as "samples.tsv".
3. Copy information about analytes.
- Click on study's tab "Show named metabolites".
- Save file as "analytes.tsv".
4. Copy information about measurements.
- Click on study's tab "Download named metabolite data".
- Save file as "measurements.tsv".
5. Copy information about signals.
- Click on study's tab "Download all metabolite data".
- Save file as "signals.tsv".
6. Correct errors.
- Project PR000305, Study ST000390
-- _"analytes.tsv"_
--- _"ascorbic acid"_, _"Vitamin C"_, change PubChem identifier from _"5785"_
to _"54670067"_.
- Project PR000322, Study ST000412
-- _"analytes.tsv"_
--- _"putrescine"_, _"Putrescine"_, change PubChem identifier from _"1049"_ to
_"1045"_.
Remove record for _"maltose"_, _"Maltose"_, PubChem _"439186"_.
That record is redundant with the record for _"cellobiose"_, PubChem _"6255"_.



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

[MetaNetX][12] is a [repository][13] of metabolic models with tools to curate
these models.

1. Import the following file to [MetaNetX][13].
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

$ metabonet clean -d /dock/

[1]: https://github.com/tcameronwaller/metabonet
[2]: https://github.com/
[3]: https://www.ncbi.nlm.nih.gov/pubmed/29078384
[4]: https://zenodo.org/record/583326
[5]: https://zenodo.org/
[6]: https://www.ncbi.nlm.nih.gov/pubmed/29140435
[7]: http://www.hmdb.ca/
[8]: https://www.ncbi.nlm.nih.gov/pubmed/26467476
[9]: http://www.metabolomicsworkbench.org/
[10]: https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000305
[11]: https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000058
[12]: https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000322
[13]: https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000599


[10]: https://www.ncbi.nlm.nih.gov/pubmed/26527720
[11]: https://www.metanetx.org/
