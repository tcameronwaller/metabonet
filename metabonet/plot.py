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
        weight=1000,
        size=25
    )
    font_three = matplotlib.font_manager.FontProperties(
        family="sans-serif",
        style="normal",
        variant="normal",
        stretch=500,
        weight=1000,
        size=15
    )
    # Compile and return references.
    return {
        "font_one": font_one,
        "font_two": font_two,
        "font_three": font_three
    }


def plot_degrees(
    nodes_metabolites_completion=None,
    nodes_metabolites_simplification=None,
    fonts=None
):
    """
    Creates a chart to represent distributions of degrees of metabolites'
    nodes.

    arguments:
        nodes_metabolites_completion (dict<dict>): information about
            metabolites' nodes
        nodes_metabolites_simplification (dict<dict>): information about
            metabolites' nodes
        fonts (dict<object>): references to definitions of font properties

    raises:

    returns:
        (object): chart object

    """

    series_one = utility.collect_value_from_records(
        key="degree",
        records=nodes_metabolites_completion.values()
    )
    series_two = utility.collect_value_from_records(
        key="degree",
        records=nodes_metabolites_simplification.values()
    )
    chart_degrees = plot_two_distributions_histograms(
        series_one=series_one,
        series_two=series_two,
        name_one="completion",
        name_two="simplification",
        fonts=fonts
    )
    # Return reference to figure.
    return chart_degrees


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


# TODO: this old version might be useful for word clouds?

def prepare_metabolites_ranks_maybe_useful_for_word_clouds(
    nodes_metabolites=None
):
    """
    Prepares summary information about ranks of metabolites' nodes.

    arguments:
        nodes_metabolites (dict<dict>): information about metabolites' nodes

    raises:

    returns:
        (dict<list<dict<str>>): summaries of metabolites' nodes' ranks

    """

    # TODO: do not need to split the ranks up... keep together...

    # Filter records to include the nodes with the top 25 total ranks.
    ranks_total = extract_sort_filter_ranks(
        rank="rank",
        filter_by="count",
        count=26,
        nodes_metabolites=nodes_metabolites
    )
    # Filter records to include all nodes with the top 25 total ranks.
    ranks_total_identifiers = utility.collect_value_from_records(
        key="identifier",
        records=ranks_total
    )
    ranks_degree = extract_sort_filter_ranks(
        rank="rank_centrality_degree",
        filter_by="identifier",
        identifiers=ranks_total_identifiers,
        nodes_metabolites=nodes_metabolites
    )
    ranks_betweenness = extract_sort_filter_ranks(
        rank="rank_centrality_betweenness",
        filter_by="identifier",
        identifiers=ranks_total_identifiers,
        nodes_metabolites=nodes_metabolites
    )
    # Compile and return information.
    return {
        "total": ranks_total,
        "degree": ranks_degree,
        "betweenness": ranks_betweenness
    }


def extract_sort_filter_ranks(
    rank=None,
    filter_by=None,
    count=None,
    identifiers=None,
    nodes_metabolites=None
):
    """
    Extracts, sorts, and filters information about ranks of metabolites' nodes.

    arguments:
        rank (str): key of rank's value in each entry
        filter_by (str): whether to filter records by count of ranks or by
            identifiers
        count (int): count of ranks to include in summary
        identifiers (list<str>): identifiers of nodes to include in summary
        nodes_metabolites (dict<dict>): information about metabolites' nodes

    raises:

    returns:
        (list<dict<str>): summary of metabolites' nodes' ranks

    """

    # Extract relevant information.
    ranks_raw = []
    for node in nodes_metabolites.values():
        identifier = node["identifier"]
        name = node["name"]
        rank_node = node[rank]
        record = {
            "identifier": identifier,
            "name": name,
            "rank": rank_node
        }
        ranks_raw.append(record)
    # Sort information by rank.
    ranks_sort = sorted(
        ranks_raw,
        key=lambda record: record["rank"],
        reverse=False
    )
    # Filter information.
    if filter_by == "count":
        ranks_filter = list(filter(
            lambda record: record["rank"] < count,
            ranks_sort
        ))
    elif filter_by == "identifier":
        ranks_filter = list(filter(
            lambda record: record["identifier"] in identifiers,
            ranks_sort
        ))
    return ranks_filter


def determine_ranks_summary_extremes(ranks_summary=None):
    """
    Determines minimal and maximal values of all ranks in summary.

    arguments:
        ranks_summary (list<dict<str>): summary of metabolites' nodes' ranks

    raises:

    returns:
        (tuple<int, int>): minimal and maximal ranks

    """

    ranks = []
    for record in ranks_summary:
        ranks.append(record["rank_total"])
        ranks.append(record["rank_degree"])
        ranks.append(record["rank_betweenness"])
    minimum = min(ranks)
    maximum = max(ranks)
    return (minimum, maximum)


def plot_ranks(
    count=None,
    nodes_metabolites=None,
    fonts=None
):
    """
    Creates a chart to represent ranks of metabolites' nodes.

    arguments:
        count (int): count of metabolites' nodes to include in summary
        nodes_metabolites (dict<dict>): information about metabolites' nodes
        fonts (dict<object>): references to definitions of font properties

    raises:

    returns:
        (object): chart object

    """

    ranks_summary = prepare_ranks_summary(
        count=count,
        nodes_metabolites=nodes_metabolites
    )
    chart_ranks = plot_three_ranks_parallel_coordinates(
        ranks_summary=ranks_summary,
        fonts=fonts
    )
    # Return reference to figure.
    return chart_ranks


def prepare_ranks_summary(
    count=None,
    nodes_metabolites=None
):
    """
    Prepares summary information about ranks of metabolites' nodes.

    arguments:
        count (int): count of metabolites' nodes to include in summary
        nodes_metabolites (dict<dict>): information about metabolites' nodes

    raises:

    returns:
        (list<dict<str>): summary of metabolites' nodes' ranks

    """

    # Extract relevant information.
    summary_raw = []
    for node in nodes_metabolites.values():
        identifier = node["identifier"]
        name = node["name"]
        rank_total = node["rank"]
        rank_degree = node["rank_centrality_degree"]
        rank_betweenness = node["rank_centrality_betweenness"]
        record = {
            "identifier": identifier,
            "name": name,
            "rank_total": rank_total,
            "rank_degree": rank_degree,
            "rank_betweenness": rank_betweenness
        }
        summary_raw.append(record)
    # Sort summary by total rank.
    summary_sort = sorted(
        summary_raw,
        key=lambda record: record["rank_total"],
        reverse=False
    )
    # Filter summary by total rank.
    summary_filter = list(filter(
        lambda record: record["rank_total"] < (count + 1),
        summary_sort
    ))
    return summary_filter


def plot_three_ranks_parallel_coordinates(
    ranks_summary=None,
    fonts=None
):
    """
    Creates a parallel coordinates chart to represent ranks of metabolites'
    nodes.

    arguments:
        ranks_summary (list<dict<str>): summary of metabolites' nodes' ranks
        fonts (dict<object>): references to definitions of font properties

    raises:

    returns:
        (object): chart object

    """

    # Determine minimal and maximal values for all axes.
    minimum, maximum = determine_ranks_summary_extremes(
        ranks_summary=ranks_summary
    )
    # Extract information for each axis.
    categories = ["Degree Rank", "Total Rank", "Betweenness Rank"]
    # Create chart.
    figure = pyplot.figure(
        figsize=(11.811, 15.748),
        tight_layout=True
    )
    axes = pyplot.axes()
    axes.invert_yaxis()
    axes.set_ylabel(
        ylabel="Node Ranks",
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
    for category in categories:
        axes.axvline(
            category,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=(0.0, 0.0, 0.0, 1.0)
        )
    # Create a plot for each value.
    for record in ranks_summary:
        values = [
            record["rank_degree"],
            record["rank_total"],
            record["rank_betweenness"]
        ]
        axes.plot(
            categories,
            values,
            color=(0.0, 0.2, 0.5, 0.75),
            linewidth=1.0
        )
        axes.text(
            categories[1],
            values[1],
            record["name"],
            backgroundcolor=(1.0, 1.0, 1.0, 0.5),
            color=(0.0, 0.0, 0.0, 1.0),
            fontproperties=fonts["font_three"],
            horizontalalignment="center",
            verticalalignment="center"
        )
    # Return reference to figure.
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
    path_degrees = os.path.join(path, "metabolite_degrees.svg")
    path_ranks_completion = os.path.join(
        path, "metabolite_ranks_completion.svg"
    )
    path_ranks_simplification = os.path.join(
        path, "metabolite_ranks_simplification.svg"
    )
    # Write information to file.
    information["chart_degrees"].savefig(
        path_degrees,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["chart_ranks_completion"].savefig(
        path_ranks_completion,
        format="svg",
        dpi=600,
        facecolor="w",
        edgecolor="w"
    )
    information["chart_ranks_simplification"].savefig(
        path_ranks_simplification,
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
    # Plot distributions for nodes' degrees before and after simplification.
    chart_degrees = plot_degrees(
        nodes_metabolites_completion=source["nodes_metabolites_completion"],
        nodes_metabolites_simplification=(
            source["nodes_metabolites_simplification"]
        ),
        fonts=fonts
    )
    # Parallel coordinates chart for ranks of metabolites' nodes.
    chart_ranks_completion = plot_ranks(
        count=20,
        nodes_metabolites=source["nodes_metabolites_completion"],
        fonts=fonts
    )
    chart_ranks_simplification = plot_ranks(
        count=20,
        nodes_metabolites=source["nodes_metabolites_simplification"],
        fonts=fonts
    )
    # Word cloud.
    # TODO: implement...
    chart_names_completion = plot_names(
        count=500,
        nodes_metabolites=source["nodes_metabolites_completion"],
        fonts=fonts
    )
    chart_names_simplification = plot_names(
        count=500,
        nodes_metabolites=source["nodes_metabolites_simplification"],
        fonts=fonts
    )

    # Compile information.
    information = {
        "chart_degrees": chart_degrees,
        "chart_ranks_completion": chart_ranks_completion,
        "chart_ranks_simplification": chart_ranks_simplification
    }
    #Write product information to file.
    write_product(directory=directory, information=information)
