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
import json

# Relevant

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
    path = os.path.join(directory, "network")
    path_nodes_reactions = os.path.join(path, "nodes_reactions.pickle")
    path_nodes_metabolites = os.path.join(path, "nodes_metabolites.pickle")
    path_links = os.path.join(path, "links.pickle")
    # Read information from file.
    with open(path_nodes_reactions, "rb") as file_source:
        nodes_reactions = pickle.load(file_source)
    with open(path_nodes_metabolites, "rb") as file_source:
        nodes_metabolites = pickle.load(file_source)
    with open(path_links, "rb") as file_source:
        links = pickle.load(file_source)
    # Compile and return information.
    return {
        "nodes_reactions": nodes_reactions,
        "nodes_metabolites": nodes_metabolites,
        "links": links
    }


def convert_networkx(
    nodes_reactions=None,
    nodes_metabolites=None,
    links=None
):
    """
    Converts information about network's nodes and links to format for
    NetworkX.

    arguments:
        nodes_reactions (dict<dict>): information about reactions' nodes
        nodes_metabolites (dict<dict>): information about metabolites' nodes
        links (dict<dict>): information about links between nodes for reactions
            and metabolites

    raises:

    returns:
        (dict<list<tuple>>): information about network's nodes and links

    """

    nodes_networkx = convert_nodes_networkx(
        nodes_reactions=nodes_reactions,
        nodes_metabolites=nodes_metabolites
    )
    links_networkx = convert_links_networkx(links=links)
    # Compile and return information.
    return {
        "nodes": nodes_networkx,
        "links": links_networkx
    }


def convert_nodes_networkx(
    nodes_reactions=None,
    nodes_metabolites=None
):
    """
    Converts information about network's nodes to format for NetworkX.

    arguments:
        nodes_reactions (dict<dict>): information about reactions' nodes
        nodes_metabolites (dict<dict>): information about metabolites' nodes

    raises:

    returns:
        (list<tuple>): information about network's nodes

    """

    nodes_networkx = []
    for node_reaction in nodes_reactions.values():
        node_networkx = (node_reaction["identifier"], node_reaction)
        nodes_networkx.append(node_networkx)
    for node_metabolite in nodes_metabolites.values():
        node_networkx = (node_metabolite["identifier"], node_metabolite)
        nodes_networkx.append(node_networkx)
    return nodes_networkx


def convert_links_networkx(links=None):
    """
    Converts information about network's links to format for NetworkX.

    arguments:
        links (dict<dict>): information about links between nodes for reactions
            and metabolites

    raises:

    returns:
        (list<tuple>): information about network's links

    """

    links_networkx = []
    for link in links.values():
        link_networkx = (link["source"], link["target"], link)
        links_networkx.append(link_networkx)
    return links_networkx


def convert_cytoscape(
    nodes_reactions=None,
    nodes_metabolites=None,
    links=None
):
    """
    Converts information about network's nodes and links to format for
    CytoScape.

    arguments:
        nodes_reactions (dict<dict>): information about reactions' nodes
        nodes_metabolites (dict<dict>): information about metabolites' nodes
        links (dict<dict>): information about links between nodes for reactions
            and metabolites

    raises:

    returns:
        (dict<dict<list<dict>>>): information about network's nodes and links

    """

    nodes_cytoscape = convert_nodes_cytoscape(
        nodes_reactions=nodes_reactions,
        nodes_metabolites=nodes_metabolites
    )
    links_cytoscape = convert_links_cytoscape(links=links)
    # Compile and return information.
    return {
        "elements": {
            "nodes": nodes_cytoscape,
            "edges": links_cytoscape
        }
    }


def convert_nodes_cytoscape(
    nodes_reactions=None,
    nodes_metabolites=None
):
    """
    Converts information about network's nodes to format for CytoScape.

    arguments:
        nodes_reactions (dict<dict>): information about reactions' nodes
        nodes_metabolites (dict<dict>): information about metabolites' nodes

    raises:

    returns:
        (list<dict>): information about network's nodes

    """

    nodes_cytoscape = []
    for node_reaction in nodes_reactions.values():
        node_reaction["id"] = node_reaction["identifier"]
        node_cytoscape = {
            "data": node_reaction
        }
        nodes_cytoscape.append(node_cytoscape)
    for node_metabolite in nodes_metabolites.values():
        node_metabolite["id"] = node_metabolite["identifier"]
        node_cytoscape = {
            "data": node_metabolite
        }
        nodes_cytoscape.append(node_cytoscape)
    return nodes_cytoscape


def convert_links_cytoscape(links=None):
    """
    Converts information about network's links to format for CytoScape.

    arguments:
        links (dict<dict>): information about links between nodes for reactions
            and metabolites

    raises:

    returns:
        (list<dict>): information about network's links

    """

    links_cytoscape = []
    for link in links.values():
        link["id"] = link["identifier"]
        link_cytoscape = {
            "data": link
        }
        links_cytoscape.append(link_cytoscape)
    return links_cytoscape


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
    path = os.path.join(directory, "conversion")
    utility.confirm_path_directory(path)
    path_networkx = os.path.join(path, "networkx.pickle")
    path_cytoscape = os.path.join(path, "cytoscape.json")
    # Write information to file.
    with open(path_networkx, "wb") as file_product:
        pickle.dump(information["networkx"], file_product)
    with open(path_cytoscape, "w") as file_product:
        json.dump(information["cytoscape"], file_product)


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to convert information about network's
    elements to versatile formats, specifically for compatibility with
    NetworkX and CytoScape.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # Convert information format for export to NetworkX.
    networkx = convert_networkx(
        nodes_reactions=source["nodes_reactions"],
        nodes_metabolites=source["nodes_metabolites"],
        links=source["links"]
    )
    # Convert information format for export to CytoScape.
    cytoscape = convert_cytoscape(
        nodes_reactions=source["nodes_reactions"],
        nodes_metabolites=source["nodes_metabolites"],
        links=source["links"]
    )
    # Compile information.
    information = {
        "networkx": networkx,
        "cytoscape": cytoscape
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
