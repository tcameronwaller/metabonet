"""
Remove intermediate directories and files from metabocurator's procedures.

Title:

    clean

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
    argparse: Package to interpret user parameters from terminal.
    textwrap: Package to format text.
    csv: Package to organize information in text.
    copy: Package to copy objects.
    pickle: Package to preserve information.
    numpy: Package to calculate with arrays of numbers.
    pandas: Package to organize collections of variables.

Classes:

    This module does not contain any classes.

Exceptions:

    This module does not contain any exceptions.

Functions:

    ...

Author:

    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    University of Utah
    Room 5520C, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

License:

    This file is part of project metabonet
    (https://github.com/tcameronwaller/metabonet/).

    MetaboNet supports custom definition of metabolic networks.
    Copyright (C) 2018 Thomas Cameron Waller

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.
    If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard
import os

# Relevant

# Custom

#dir()
#importlib.reload()

###############################################################################
# Functionality


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to convert information about metabolic
    entities and sets to versatile formats.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Specify directories and files.
    path_reconciliation = os.path.join(directory, "reconciliation")
    path_extraction = os.path.join(directory, "extraction")
    path_extrication = os.path.join(directory, "extrication")
    path_enhancement = os.path.join(directory, "enhancement")
    path_curation = os.path.join(directory, "curation")
    # Remove directories and files.
    if (False):
        os.remove(path_reconciliation, "recon2m2_reconciliation.xml")
        os.remove(path_reconciliation, "recon2m2_metanetx_compartments.tsv")
        os.remove(path_reconciliation, "recon2m2_metanetx_genes.tsv")
        os.remove(path_reconciliation, "recon2m2_metanetx_metabolites.tsv")
        os.remove(path_reconciliation, "recon2m2_metanetx_reactions.tsv")
    os.remove(os.path.join(path_extraction, "compartments.pickle"))
    os.remove(os.path.join(path_extraction, "processes.pickle"))
    os.remove(os.path.join(path_extraction, "reactions.pickle"))
    os.remove(os.path.join(path_extraction, "metabolites.pickle"))
    os.remove(os.path.join(path_extraction, "reactions.tsv"))
    os.remove(os.path.join(path_extraction, "metabolites.tsv"))
    os.rmdir(path_extraction)
    os.remove(os.path.join(path_extrication, "hmdb_summary.pickle"))
    os.remove(os.path.join(path_extrication, "hmdb_summary.tsv"))
    os.rmdir(path_extrication)
    os.remove(os.path.join(path_enhancement, "compartments.pickle"))
    os.remove(os.path.join(path_enhancement, "processes.pickle"))
    os.remove(os.path.join(path_enhancement, "reactions.pickle"))
    os.remove(os.path.join(path_enhancement, "metabolites.pickle"))
    os.remove(os.path.join(path_enhancement, "reactions.tsv"))
    os.remove(os.path.join(path_enhancement, "metabolites.tsv"))
    os.remove(os.path.join(path_enhancement, "reactions_filter.tsv"))
    os.rmdir(path_enhancement)
    os.remove(os.path.join(path_curation, "compartments.pickle"))
    os.remove(os.path.join(path_curation, "processes.pickle"))
    os.remove(os.path.join(path_curation, "reactions.pickle"))
    os.remove(os.path.join(path_curation, "metabolites.pickle"))
    os.remove(os.path.join(path_curation, "reactions.tsv"))
    os.remove(os.path.join(path_curation, "metabolites.tsv"))
    os.rmdir(path_curation)
