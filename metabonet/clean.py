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
import utility

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
    path_candidacy = os.path.join(directory, "candidacy")
    path_network = os.path.join(directory, "network")
    path_conversion = os.path.join(directory, "conversion")
    path_measurement = os.path.join(directory, "measurement")
    path_analysis = os.path.join(directory, "analysis")
    path_plot = os.path.join(directory, "plot")
    # Remove directories and files.
    # Candidacy.
    utility.remove_file(os.path.join(path_candidacy, "metabolites.pickle"))
    utility.remove_file(os.path.join(path_candidacy, "reactions.pickle"))
    utility.remove_file(os.path.join(path_candidacy, "reactions.tsv"))
    utility.remove_file(os.path.join(path_candidacy, "metabolites.tsv"))
    utility.remove_empty_directory(path_candidacy)
    # Network.
    utility.remove_file(os.path.join(path_network, "links.pickle"))
    utility.remove_file(os.path.join(path_network, "nodes_metabolites.pickle"))
    utility.remove_file(os.path.join(path_network, "nodes_reactions.pickle"))
    utility.remove_empty_directory(path_network)
    # Conversion.
    utility.remove_file(os.path.join(
        path_conversion, "network_elements_cytoscape.json"
    ))
    utility.remove_file(os.path.join(
        path_conversion, "network_elements_networkx.pickle"
    ))
    utility.remove_empty_directory(path_conversion)
    # Measurement.
    utility.remove_file(os.path.join(path_measurement, "metabolites.pickle"))
    utility.remove_file(os.path.join(path_measurement, "reactions.pickle"))
    utility.remove_file(os.path.join(path_measurement, "reactions.tsv"))
    utility.remove_file(os.path.join(path_measurement, "metabolites.tsv"))
    utility.remove_file(os.path.join(path_measurement, "metabolites.tsv"))
    utility.remove_empty_directory(path_measurement)
    # Analysis.
    utility.remove_file(os.path.join(path_analysis, "nodes_reactions.tsv"))
    utility.remove_file(os.path.join(path_analysis, "nodes_metabolites.tsv"))
    utility.remove_file(os.path.join(path_analysis, "network_reactions.tsv"))
    utility.remove_file(os.path.join(path_analysis, "network_metabolites.tsv"))
    utility.remove_empty_directory(path_analysis)
    # Plot.
    utility.remove_file(os.path.join(path_plot, "measurements_one.svg"))
    utility.remove_file(os.path.join(path_plot, "measurements_two.svg"))
    utility.remove_file(os.path.join(path_plot, "measurements_three.svg"))
    utility.remove_file(os.path.join(path_plot, "measurements_four.svg"))
    utility.remove_file(os.path.join(path_plot, "measurements_five.svg"))
    utility.remove_empty_directory(path_plot)
