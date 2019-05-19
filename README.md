# MetaboNet

This project supports curation of the human metabolic models, definition of custom metabolic networks, and analyses of these networks.

## Launch MetaboNet's Singularity Container

Coming soon...
```bash
$ docker pull jordanberg/metabonet:publication
$ docker run metabonet <insert commands here>
```

## Install MetaboNet

Download [MetaboNet][1] from [GitHub][3].

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

# View MetaboNet sub-module help menus

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
MetaboNet repository installation.

## Set Up Source Information Within Dock

## Metabolic Information

MetaboNet derives information from both a model of human metabolism and from a
metabolite database. When installing MetaboNet, these files are automatically installed;
however, if you wish to use different version, you can follow the instructions below, substituting
file names for the appropriate versioned file. As the necessary files are large,
these instructions specify how to access them from original repositories.

### Download Model of Human Metabolism

MetaboNet derives information about human metabolism from model
[Recon 2M.2](4). Access [Recon 2M.2](4) by downloading file
`Recon2M.2_MNX_Entrez_Gene.xml` (14.2 MB) from [Zenodo](6) record [583326](5).
Follow the directions below to change the name of this file to `recon2m2.xml`
and place it in the following location.

```bash
$ cd /path/to/metabonet/
$ curl -O https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml?download=1
$ mv Recon2M.2_MNX_Entrez_Gene.xml?download=1 dock/source/recon2m2.xml
```

### Download Human Metabolite Database

MetaboNet enhances information about human metabolites from the
[Human Metabolome Database](7) [HMDB](8). Access the [HMDB](8) by
downloading version `4.0` of file `hmdb_metabolites.xml` (4.2 GB).
Decompress this file and place it in the following location, as shown below.

```bash
$ cd /path/to/metabonet/
$ curl -O http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
$ unzip hmdb_metabolites.zip
$ mv hmdb_metabolites.xml dock/source/
$ rm hmdb_metabolites.zip
```

## Customization Parameters

MetaboNet is versatile in its support for customization at many levels. Files
of default customization parameters are available within the "customization"
directory in the [MetaboNet](1) repository. These are already provided, but can
removed by doing the following:

```bash
$ cd /path/to/metabonet/
$ rm dock/source/customization/*
```

If you wish to edit these parameters, please open the appropriate file and modify the relevant
parameters with your favorite text editor.

## Metabolomic Measurements

MetaboNet curates and analyzes metabolomic measurements. The original data sets
come from [Metabolomics Workbench](10), a [repository](11) for data from metabolomic
studies. For the example case, the information is already provided; however, if one wishes
to use a custom metabolomic dataset, the user should follow the steps below.

1. Locate specific project and study on [Metabolomics Workbench][11].
2. Copy information about samples.
  - Click on study's tab `Show all samples`.
  - Copy information from sample table.
  - Save file as `samples.tsv` to `dock/source/measurements/`.
3. Copy information about analytes.
  - Click on study's tab `Show named metabolites`.
  - Save file as `analytes.tsv`.
4. Copy information about measurements.
  - Click on study's tab `Download named metabolite data`.
  - Save file as "measurements.tsv".
5. Copy information about signals.
  - Click on study's tab `Download all metabolite data`.
  - Save file as `signals.tsv`.
6. Correct errors.
  - `Project PR000305/Study ST000390/analytes.tsv`
    - `"ascorbic acid", "Vitamin C"`, change PubChem identifier from `5785` to `54670067`.
  - `Project PR000322/Study ST000412/analytes.tsv`
    - `"putrescine", "Putrescine"`, change PubChem identifier from `1049` to `1045`.
    - Remove record for `"maltose", "Maltose", PubChem 439186` (Record is redundant with the record for `"cellobiose", PubChem 6255`)


## Curate Metabolic Model and Metabolomic Measurements

MetaboNet's `model` routine includes multiple procedures to curate the
metabolic model and to curate metabolomic measurements. This curation prepares
to define custom metabolic networks and integrate measurements with these
networks. NOTE: The next two steps are not required if using the default Recon and HMDB versions
provided during MetaboNet installation.

### (Optional) Reconcile Metabolic Model for Import to MetaNetX

Execute the `reconciliation` procedure of the `model` routine in MetaboNet.

```bash
$ metabonet model -d dock/ -r
```

This command should take a couple seconds to execute, after which you should see output similar to below:

```bash
--------------------------------------------------
... call to model routine ...
... executing reconciliation procedure ...
compartments: 10
reactions: 5842
metabolites: 4000
```

### (Optional) Retrieve Metabolic Model Files from MetaNetX
[MetaNetX](13) is a repository of metabolic models with tools to curate
these models.

1. Import the following file to [MetaNetX](https://www.metanetx.org/cgi-bin/mnxweb/import_mnet).

```bash
dock/reconciliation/recon2m2_reconciliation.xml
```

2. Once the import has finished (may take a few minutes), click on the hyperlink in the import table named after the name you provided before importing. This will open a table view of the import results. Scroll down and download the following files within the `Model parts` row to the `dock/reconciliation` directory with the MetaboNet directory and and change the names as specified below.

```bash
compartments.tsv -> dock/reconciliation/recon2m2_metanetx_compartments.tsv
enzymes.tsv -> dock/reconciliation/recon2m2_metanetx_genes.tsv
chemicals.tsv -> dock/reconciliation/recon2m2_metanetx_metabolites.tsv
reactions.tsv -> dock/reconciliation/recon2m2_metanetx_reactions.tsv
```

3. Save the "Mapping summary". This summary can be found by clicking the `Import` hyperlink on the `Analysis` row within the table view of the import results. Then right-click the `Mapping summary` hyperlink and save the file as specified below:

```bash
dock/reconciliation/metanetx_import_report.tsv
```

This summary is useful for review and also includes useful information to curate names of metabolites.
Import category "chemical compounds mapped to the MNXref namespace with
different ID and different description" is of particular interest for curation
of information about metabolites.


### Curate the Metabolic Model

1. Execute the `collection`, `extraction`, `enhancement`, `curation`, and
`conversion` procedures of the `model` routine in MetaboNet using the command provided below:

```bash
$ metabonet model -d /dock/ -ceauv
```

The `enhancement` procedure alone requires about 1 hour to complete. The process altogether may take a couple of hours to complete. Below are example system specs:

```
System:
$ /usr/sbin/system_profiler SPHardwareDataType
Hardware:

    Hardware Overview:

      Model Name: MacBook Pro
      Model Identifier: MacBookPro15,1
      Processor Name: Intel Core i9
      Processor Speed: 2.9 GHz
      Number of Processors: 1
      Total Number of Cores: 6
      L2 Cache (per Core): 256 KB
      L3 Cache: 12 MB
      Hyper-Threading Technology: Enabled
      Memory: 32 GB
```
```
$ sysctl -a | grep cpu
hw.ncpu: 12
```

Below are execution specs for each sub-module of the above command:

```
Collection:
real	0m5.039s
user	0m5.238s
sys	0m0.249s
```
```
Extraction:
real	2m54.355s
user	2m52.744s
sys	0m1.683s
```
```
Enhancement:

```
```
Curation:

```
```
Conversion:

```



2. (Optional) Edit file of special use for automatic curation of reactions:

```bash
dock/enhancement/reactions_filter.tsv
```

3. (Optional) Customize curation by editing the following files of parameters.

```bash
dock/source/customization/curation_compartments.tsv
dock/source/customization/curation_processes.tsv
dock/source/customization/curation_reactions.tsv
dock/source/customization/curation_metabolites.tsv
```

### Curate Metabolomic Measurements

1. Execute the `measurement` procedure of the `model` routine in MetaboNet.

```bash
$ metabonet model -d dock/ -m
```

Below are execution specs for the above `measurement` command:
```

```

2. Review the following files to check the matches between analytes and
metabolites. Discrepancies might justify modifications to curation of
metabolites.

```bash
root/measurement/study_one_report.tsv
root/measurement/study_two_report.tsv
root/measurement/study_three_report.tsv
root/measurement/study_four_report.tsv
root/measurement/study_five_report.tsv
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

Execute the `candidacy`, `network`, and `conversion` procedures of the
`network` routine in MetaboNet.

- Set `compartments == True` and `hubs == True`:
```bash
$ metabonet network -d dock/ -yc -np -v
```

- Set `compartments == True` and `hubs == False`:
```bash
$ metabonet network -d /dock/ -ycs -np -v
```

- Set `compartments == False` and `hubs == True`:
```bash
$ metabonet network -d /dock/ -y -np -v
```

- Set `compartments == False` and `hubs == False`:
```bash
$ metabonet network -d /dock/ -ys -np -v
```

Below are execution specs for the `compartments == True` and `hubs == True` execution:
```

```

## Analyze Metabolic Networks

Execute the `analysis` procedure of the `network` routine in MetaboNet.

```bash
$ metabonet network -d dock/ -a
```

Run time of `analysis` procedure can be 0.5-2.0 hours depending on the network.

Below are execution specs for the `compartments == True` and `hubs == True` execution:
```

```

## Visualize Metabolic Networks

The resulting network database can be found from the following file:
```bash
dock/conversion/network_elements_cytoscape.json
```

## Integrate Metabolomic Measurements in Metabolic Network

Execute the `measurement` procedure of the `network` routine in MetaboNet.

```bash
$ metabonet network -d dock/ -m
```

The following file includes measurements from all studies ready for integration
with the metabolic network in Cytoscape.

```bash
root/measurement/metabolites.tsv
```

## Remove of automatically-generate files and directories

```bash
$ metabonet clean -d dock/
```





[1]: (https://github.com/tcameronwaller/metabonet)
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
