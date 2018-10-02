# metabonet

This project supports procedures to curate metabolic models and, from these, to
define metabolic networks of reactions between metabolites.

## Description

This document provides a description of general procedures. Information in more
detail about functionality is accessible in the help for the interface.

## Interface

The current interface for the package is a module that behaves like a script.
It is necessary to execute this file as a python script. The interface includes
multiple routines with their own parameters and help information.

$ python3 interface.py -help
$ python3 interface.py model -help
$ python3 interface.py network -help
$ python3 interface.py clean -help

## Source information

All procedures require specification of a path to a root directory, and this
directory will usually be the same for all procedures. Procedures read source
information from directories and files within this root directory. Procedures
also write product information to directories and files within this root
directory.

$ python3 interface.py {model,network,clean} -d root/directory

## Curation of metabolic model and metabolomic measurements

This procedure involves the "model" routine within MetaboNet.

$ python3 interface.py model -help
$ python3 interface.py model -d root/directory -rceauvm

1. Reconcile raw metabolic model for integration to MetaNetX.

$ python3 interface.py model -reconciliation

2. Collect information about compartments, processes, reactions, and
metabolites.

$ python3 interface.py model -collection

3. Extract reference information about metabolites from the Human Metabolome
Database (HMDB).

$ python3 interface.py model -extraction

4. Enhance information about metabolites and reactions.

$ python3 interface.py model -enhancement

5. Curate information about individual compartments, processes, reactions, and
metabolites.

$ python3 interface.py model -curation

6. Convert information for review and export.

$ python3 interface.py model -conversion

7. Curate metabolomic measurements from multiple studies.

$ python3 interface.py model -measurement

## Definition of custom metabolic networks and association of measurements

$ python3 interface.py network -help


## Removal of automatically-generate files and directories

$ python3 interface.py clean -help
