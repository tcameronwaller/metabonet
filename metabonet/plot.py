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

# Relevant.
import numpy
import matplotlib
matplotlib.use('agg')
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


def define_font_properties():
    """
    Defines font properties.

    arguments:

    raises:

    returns:
        (dict<object>): references to definitions of font properties

    """

    font_one = matplotlib.font_manager.FontProperties(
        family="sans-serif",
        style="normal",
        variant="normal",
        stretch=500,
        weight=1000,
        size=30
    )
    font_two = matplotlib.font_manager.FontProperties(
        family="sans-serif",
        style="normal",
        variant="normal",
        stretch=500,
        weight=500,
        size=25
    )
    # Compile and return references.
    return {
        "font_one": font_one,
        "font_two": font_two,
    }


def plot_two_distributions_histograms(
    series_one=None,
    series_two=None,
    name_one=None,
    name_two=None,
    fonts=None
):
    """
    Creates a histogram chart to represent frequency distributions of two
    series.

    arguments:
        series_one (list<int>): series of counts
        series_two (list<int>): series of counts
        name_one (str): name of distribution
        name_two (str): name of distribution
        fonts (dict<object>): references to definitions of font properties

    raises:

    returns:
        (object): chart object

    """

    # Determine histogram bins.
    hist, bin_edges = numpy.histogram(series_one, bins=50)
    # Create chart.
    figure = pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = pyplot.axes()
    values_one, bins_one, patches_one = axes.hist(
        series_one,
        bins=bin_edges,
        histtype="bar",
        align="left",
        orientation="vertical",
        rwidth=0.35,
        log=True,
        color=(0.0, 0.2, 0.5, 1.0),
        label=name_one,
        stacked=False
    )
    values_two, bins_two, patches_two = axes.hist(
        series_two,
        bins=bin_edges,
        histtype="bar",
        align="mid",
        orientation="vertical",
        rwidth=0.35,
        log=True,
        color=(1.0, 0.6, 0.2, 1.0),
        label=name_two,
        stacked=False
    )
    axes.legend(
        loc="upper right",
        markerscale=2.5,
        markerfirst=True,
        prop=fonts["font_one"],
        edgecolor=(0.0, 0.0, 0.0, 1.0)
    )
    axes.set_xlabel(
        xlabel="Node Degree (Count)",
        labelpad=25,
        alpha=1.0,
        backgroundcolor=(1.0, 1.0, 1.0, 1.0),
        color=(0.0, 0.0, 0.0, 1.0),
        fontproperties=fonts["font_one"]
    )
    axes.set_ylabel(
        ylabel="Count of Nodes",
        labelpad=25,
        alpha=1.0,
        backgroundcolor=(1.0, 1.0, 1.0, 1.0),
        color=(0.0, 0.0, 0.0, 1.0),
        fontproperties=fonts["font_one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=(0.0, 0.0, 0.0, 1.0),
        pad=5,
        labelsize=fonts["font_two"].get_size(),
        labelcolor=(0.0, 0.0, 0.0, 1.0)
    )
    return figure


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
    path_histogram = os.path.join(path, "metabolite_degree_histogram.svg")
    # Write information to file.
    information["chart_histogram"].savefig(
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
    # Define font properties.
    fonts = define_font_properties()
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
        name_two="simplification",
        fonts=fonts
    )
    # Parallel coordinates plot.
    # http://benalexkeen.com/parallel-coordinates-in-matplotlib/

    # Word cloud.

    # Compile information.
    information = {
        "chart_histogram": chart_histogram,
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
