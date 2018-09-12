"""
Evaluate candidacy of metabolites and reactions for representation in network.

Candidate entities, metabolites and reactions, are entities that are elligible
candidates for representation in the network.
An entity's candidacy depends on the context of interest.
This context includes relevance compartmentalization in general, relevance of
specific compartments and processes, and relevance of individual entities.
An entity's candidacy also depends on the candidacies of other entities to
which the entity relates.
A reaction's candidacy depends on the candidacies of metabolites that
participate in it.
A metabolite's candidacy depends on the candidacies of reactions in which
it participates.

Title:

    candidacy

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
    argparse: Package to interpret user parameters from terminal.
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
import pickle
import copy

# Relevant
import networkx as ntx

# Custom
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source(directory=None):
    """
    Reads and organizes source information from file

    arguments:
        directory (str): directory of source files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_conversion = os.path.join(directory, "conversion")
    path_compartments = os.path.join(path_conversion, "compartments.pickle")
    path_processes = os.path.join(path_conversion, "processes.pickle")
    path_reactions = os.path.join(path_conversion, "reactions.pickle")
    path_metabolites = os.path.join(path_conversion, "metabolites.pickle")
    path_candidacy = os.path.join(directory, "candidacy")
    path_reactions_candidacy = os.path.join(path_candidacy, "reactions.pickle")
    path_metabolites_candidacy = os.path.join(path_candidacy, "metabolites.pickle")
    # Read information from file.
    with open(path_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    with open(path_reactions_candidacy, "rb") as file_source:
        reactions_candidacy = pickle.load(file_source)
    with open(path_metabolites_candidacy, "rb") as file_source:
        metabolites_candidacy = pickle.load(file_source)
    # Compile and return information.
    return {
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites,
        "reactions_candidacy": reactions_candidacy,
        "metabolites_candidacy": metabolites_candidacy,
    }


def write_product(directory=None, information=None):
    """
    Writes product information to file

    arguments:
        directory (str): directory for product files
        information (object): information to write to file

    raises:

    returns:

    """

    # Specify directories and files.
    path = os.path.join(directory, "network")
    utility.confirm_path_directory(path)
    path_nodes_reactions = os.path.join(path, "nodes_reactions.pickle")
    path_nodes_metabolites = os.path.join(path, "nodes_metabolites.pickle")
    path_links = os.path.join(path, "links.pickle")
    # Write information to file.
    with open(path_nodes_reactions, "wb") as file_product:
        pickle.dump(information["nodes_reactions"], file_product)
    with open(path_nodes_metabolites, "wb") as file_product:
        pickle.dump(information["nodes_metabolites"], file_product)
    with open(path_links, "wb") as file_product:
        pickle.dump(information["links"], file_product)


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to define network's elements.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)


    network = ntx.DiGraph()
    # Method add_nodes_from accepts a list of tuples of node identifier (str),
    # node attributes (dict).
    network.add_nodes_from(nodes_records)
    # Method add_edges_from accepts a list of tuples of source node identifier
    # (str), target node identifier (str), link attributes (dict).
    network.add_edges_from(links_records)

    # list(network.nodes)
    # list(network.edges)
    # list(network.neighbors([node_1, node_2]))
    # list(network.degree([node_1, node_2]))
    # network.nodes.data()
    # list(ntx.connected_components(network))
    # ntx.clustering(network)


    if False:
        # Compile information.
        information = {
            "nodes_reactions": nodes_reactions,
            "nodes_metabolites": nodes_metabolites,
            "links": links
        }
        #Write product information to file.
        write_product(directory=directory, information=information)
