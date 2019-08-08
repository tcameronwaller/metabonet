# [MetaboNet][1]

This project supports curation of the model of human cellular metabolism,
definition of custom metabolic networks to represent this model, and analyses
of these networks.



## License and citation

This file is part of project [MetaboNet][1].

> [MetaboNet][1] supports custom definition of metabolic networks.
Copyright (C) 2018 Thomas Cameron Waller

> This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

> This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

> You should have received a copy of the GNU General Public License along
with this program.
If not, see <http://www.gnu.org/licenses/>.

If you use [MetaboNet][1], please cite its publication.

> T. Cameron Waller, Jordan A. Berg, Brian E. Chapman, Jared P. Rutter.
> "Compartments and Hubs Differentiate the Relevance of Metabolic Networks to
Metabolomic Experiments".
> Journal. Year. Issue.
> link



## Archives

Many users might find it more convenient to access archive versions of major
exports from [MetaboNet][1].

| export from [MetaboNet][1]                    | file in archive[ref] on [Zenodo][5]      | explanation                                     |
| :-------------------------------------------- | :--------------------------------------- | :---------------------------------------------- |
| ~/dock/conversion/dymetabonet.json            | ~/model/dymetabonet.json                 | metabolic model format for [DyMetaboNet][20]    |
| ~/dock/conversion/compartments.pickle         | ~/model/compartments.pickle              | information about cellular compartments         |
| ~/dock/conversion/compartments.tsv            | ~/model/compartments.tsv                 | text abbreviation of cellular compartments      |
| ~/dock/conversion/processes.pickle            | ~/model/processes.pickle                 | information about metabolic processes           |
| ~/dock/conversion/processes.tsv               | ~/model/processes.tsv                    | text abbreviation of metabolic processes        |
| ~/dock/conversion/reactions.pickle            | ~/model/reactions.pickle                 | information about chemical reactions            |
| ~/dock/conversion/reactions.tsv               | ~/model/reactions.tsv                    | text abbreviation of chemical reactions         |
| ~/dock/conversion/metabolites.pickle          | ~/model/metabolites.pickle               | information about metabolites                   |
| ~/dock/conversion/metabolites.tsv             | ~/model/metabolites.tsv                  | text abbreviation of metabolites                |
| ~/dock/measurement/metabolites.tsv            | ~/measurement/metabolites.tsv            | integration of measurements with metabolites    |
| ~/dock/network/compartments-true_hubs-true/   | ~/network/compartments-true_hubs-true/   | files for compartmental network with hubs       |
| ~/dock/network/compartments-true_hubs-false/  | ~/network/compartments-true_hubs-false/  | files for compartmental network without hubs    |
| ~/dock/network/compartments-false_hubs-true/  | ~/network/compartments-false_hubs-true/  | files for noncompartmental network with hubs    |
| ~/dock/network/compartments-false_hubs-false/ | ~/network/compartments-false_hubs-false/ | files for noncompartmental network without hubs |
| ./links.pickle                                | ./links.pickle                           | network's links                                 |
| ./nodes_reactions.pickle                      | ./nodes_reactions.pickle                 | network's nodes for reactions                   |
| ./nodes_metabolites.pickle                    | ./nodes_metabolites.pickle               | network's nodes for metabolites                 |
| ./network_cytoscape.json                      | ./network_cytoscape.json                 | network format for Cytoscape                    |
| ./network_networkx.pickle                     | ./network_networkx.pickle                | network format for NetworkX                     |



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
dock_template LICENSE metabonet README.md setup.py
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
$ cp ~/Downloads/metabonet-master/dock_template/ ~/dock/
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

#### Option 1: Access previous measurement files in _"dock"_ directory

**_DO:_** Inspect files from metabolomic measurements in _"dock"_ directory.

```bash
$ ls ~/dock/source/measurement/
```

#### Option 2: Access and prepare metabolomic measurements from scratch

**_DO:_** Access and organize data for metabolomic measurements.

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

[MetaboNet's][1] _"model"_ routine includes multiple procedures to curate the
metabolic model and to organize metabolomic measurements. This curation
prepares the model to define custom metabolic networks and integrate
measurements with these networks.

### Reconcile Metabolic Model with [MetaNetX][14]

[MetaNetX][14] is a [repository][15] of metabolic models with useful
information about metabolites and reactions.

#### Option 1: Access previous reconciliation files in _"dock"_ directory

**_DO:_** Inspect reconciliation files in _"dock"_ directory.

```bash
$ ls ~/dock/reconciliation/
```

#### Option 2: Access and prepare reconciliation files from scratch

**_DO:_** Execute the _"reconciliation"_ procedure of the _"model"_ routine.

```bash
$ metabonet model -d ~/dock/ -r
```

**_DO:_** Integrate metabolic model with [MetaNetX][15].

1. Import the following file to [MetaNetX][15].

```bash
~/dock/reconciliation/recon2m2_reconciliation.xml
```

2. Save _"Mapping summary"_ from [MetaNetX][15] as the following file. This
summary is useful for review and also includes information to curate names of
metabolites. Import category _"chemical compounds mapped to the MNXref
namespace with different ID and different description"_ is of particular
interest for curation of information about metabolites.

```bash
/dock/reconciliation/metanetx_import_report.tsv
```

3. Save files from [MetaNetX][15] to the following directories and files.

```bash
$ cp ~/Downloads/compartments.tsv ~/dock/reconciliation/recon2m2_metanetx_compartments.tsv
$ cp ~/Downloads/enzymes.tsv ~/dock/reconciliation/recon2m2_metanetx_genes.tsv
$ cp ~/Downloads/chemicals.tsv ~/dock/reconciliation/recon2m2_metanetx_metabolites.tsv
$ cp ~/Downloads/reactions.tsv ~/dock/reconciliation/recon2m2_metanetx_reactions.tsv
```

### Curate Metabolic Model

**_DO:_** Execute the _"collection"_, _"extraction"_, _"enhancement"_,
_"curation"_, and _"conversion"_ procedures of the _"model"_ routine.

The _"enhancement"_ procedure alone requires about 1 hour to complete.

```bash
$ metabonet model -d ~/dock/ -ceauv
```

**_DO:_** Customize curation by editing the following files of parameters.

```bash
ls ~/dock/source/customization/
~/dock/source/customization/curation_compartments.tsv
~/dock/source/customization/curation_processes.tsv
~/dock/source/customization/curation_reactions_custom.tsv
~/dock/source/customization/curation_reactions.tsv
~/dock/source/customization/curation_metabolites.tsv
```

A file of special use to review automatic curation of reactions follows.

```bash
~/dock/enhancement/reactions_filter.tsv
```

### Export Metabolic Model

[MetaboNet's][1] _"model"_ routine exports information from the metabolic model
to multiple formats.

```bash
ls ~/dock/conversion/
~/dock/conversion/dymetabonet.json <- this file is compatible for import to DyMetaboNet
~/dock/conversion/compartments.pickle
~/dock/conversion/compartments.tsv
~/dock/conversion/processes.pickle
~/dock/conversion/processes.tsv
~/dock/conversion/reactions.pickle
~/dock/conversion/reactions.tsv
~/dock/conversion/metabolites.pickle
~/dock/conversion/metabolites.tsv
```



### Curate Metabolomic Measurements

**_DO:_** Execute the _"measurement"_ procedure of the _"model"_ routine.

```bash
$ metabonet model -d ~/dock/ -m
```

**_DO:_** Review the following files to check the matches between analytes and
metabolites. Discrepancies might justify modifications to curation of
metabolites.

```bash
~/dock/measurement/study_one_report.tsv
~/dock/measurement/study_two_report.tsv
~/dock/measurement/study_three_report.tsv
~/dock/measurement/study_four_report.tsv
~/dock/measurement/study_five_report.tsv
```

**_DO:_** Inspect files of summaries for metabolomic measurements in _"dock"_
directory.

```bash
~/dock/measurement/study_one.tsv
~/dock/measurement/study_two.tsv
~/dock/measurement/study_three.tsv
~/dock/measurement/study_four.tsv
~/dock/measurement/study_five.tsv
```



## Define Custom Metabolic Networks

[MetaboNet's][1] _"network"_ routine includes multiple procedures to define and
analyze customizable metabolic networks.

**_DO:_** Customize definition of metabolic networks by editing the following
files of parameters.

```bash
~/dock/source/customization/filtration_compartments.tsv
~/dock/source/customization/filtration_processes.tsv
~/dock/source/customization/simplification_reactions.tsv
~/dock/source/customization/simplification_metabolites.tsv
```

**_DO:_** Execute the _"candidacy"_, _"network"_, and _"conversion"_ procedures
of the _"network"_ routine.

compartments true, hubs true
```bash
$ metabonet network -d ~/dock/ -yc -np -v
```

compartments true, hubs false
```bash
$ metabonet network -d ~/dock/ -ycs -np -v
```

compartments false, hubs true
```bash
$ metabonet network -d ~/dock/ -y -np -v
```

compartments false, hubs false
```bash
$ metabonet network -d ~/dock/ -ys -np -v
```

### Export Metabolic Networks

[MetaboNet's][1] _"network"_ routine exports information about the metabolic
network.

```bash

ls ~/dock/network/
~/dock/network/links.pickle
~/dock/network/nodes_metabolites.pickle
~/dock/network/nodes_reactions.pickle

ls ~/dock/conversion/
~/dock/conversion/network_elements_cytoscape.json <- this file is compatible for import to Cytoscape
~/dock/conversion/network_elements_networkx.json <- this file is compatible for import to NetworkX
```



## Analyze Metabolic Networks

[MetaboNet's][1] _"network"_ routine applies functionality from the
[NetworkX][16] [application][17] to describe metabolic networks.

**_DO:_** Execute the _"analysis"_ procedure of the _"network"_ routine.

```bash
$ metabonet network -d ~/dock/ -a
```

Run time of _"analysis"_ procedure can be 0.5-2.0 hours depending on the
network.



## Visualize Metabolic Networks

The [Cytoscape][18] [application][19] is useful to visualize metabolic networks
from [MetaboNet][1].

**_DO:_** Import the following file to [Cytoscape][19].

```bash
~/dock/conversion/network_elements_cytoscape.json
```



## Integrate Metabolomic Measurements in Metabolic Network

**_DO:_** Execute the _"measurement"_ procedure of the _"network"_ routine.

```bash
$ metabonet network -d ~/dock/ -m
```

**_DO:_** Import the following file of measurements from all studies for
integration with the metabolic network in [Cytoscape][19].

```bash
~/dock/measurement/metabolites.tsv
```

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
[14]: https://www.ncbi.nlm.nih.gov/pubmed/26527720
[15]: https://www.metanetx.org/
[16]: https://www.osti.gov/biblio/960616
[17]: https://networkx.github.io/
[18]: https://www.ncbi.nlm.nih.gov/pubmed/14597658
[19]: https://cytoscape.org/
[20]: https://github.com/tcameronwaller/dymetabonet
