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
import matplotlib.pyplot as pyplot

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
    path_analysis = os.path.join(directory, "analysis")
    path_nodes_metabolites_completion = os.path.join(
        path_analysis, "nodes_metabolites_completion.pickle"
    )
    path_nodes_metabolites_simplification = os.path.join(
        path_analysis, "nodes_metabolites_simplification.pickle"
    )
    # Read information from file.
    with open(path_nodes_metabolites_completion, "rb") as file_source:
        nodes_metabolites_completion = pickle.load(file_source)
    with open(path_nodes_metabolites_simplification, "rb") as file_source:
        nodes_metabolites_simplification = pickle.load(file_source)
    # Compile and return information.
    return {
        "nodes_metabolites_completion": nodes_metabolites_completion,
        "nodes_metabolites_simplification": nodes_metabolites_simplification
    }


def plot_two_distributions_histograms(
    series_one=None,
    series_two=None,
    name_one=None,
    name_two=None
):
    """
    Creates a histogram chart to represent frequency distributions of two
    series.

    arguments:
        series_one (list<int>): series of counts
        series_two (list<int>): series of counts
        name_one (str): name of distribution
        name_two (str): name of distribution

    raises:

    returns:
        (object): chart object

    """

    #figure = plt.figure()
    #figure.suptitle("title")
    pyplot.hist(series_one, "auto", label=name_one)
    pyplot.hist(series_two, "auto", label=name_two)
    pyplot.legent(loc="upper right")




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
    path = os.path.join(directory, "plot")
    utility.confirm_path_directory(path)
    path_histogram = os.path.join(path, "metabolite_degree_histogram.pdf")
    # Write information to file.
    pyplot.savefig(
        path_histogram,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )


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
    print(source["nodes_metabolites"])

    # Histogram plot for nodes' degree distributions.
    # Plot distributions for nodes' degrees before and after simplification.
    # TODO: prepare series within a separate function?
    series_one = utility.collect_value_from_records(
        key="degree",
        records=source["nodes_metabolites_completion"].values()
    )
    series_two = utility.collect_value_from_records(
        key="degree",
        records=source["nodes_metabolites_simplification"].values()
    )
    chart_histogram = plot_two_distributions_histograms(
        series_one=series_one,
        series_two=series_two,
        name_one="completion",
        name_two="simplification"
    )
    # Parallel coordinates plot.

    # Word cloud.

    # Compile information.
    information = {
        "chart_histogram": chart_histogram,
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
