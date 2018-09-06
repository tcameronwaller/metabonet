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

# Relevant

# Custom

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
    path_customization = os.path.join(directory, "customization")
    path_filtration_compartments = os.path.join(
        path_customization, "filtration_compartments.tsv"
    )
    path_filtration_processes = os.path.join(
        path_customization, "filtration_processes.tsv"
    )
    path_simplification_reactions = os.path.join(
        path_customization, "simplification_reactions.tsv"
    )
    path_simplification_metabolites = os.path.join(
        path_customization, "simplification_metabolites.tsv"
    )
    path = os.path.join(directory, "conversion")
    path_compartments = os.path.join(path, "compartments.pickle")
    path_processes = os.path.join(path, "processes.pickle")
    path_reactions = os.path.join(path, "reactions.pickle")
    path_metabolites = os.path.join(path, "metabolites.pickle")
    # Read information from file.
    with open(path_compartments, "rb") as file_source:
        compartments = pickle.load(file_source)
    with open(path_processes, "rb") as file_source:
        processes = pickle.load(file_source)
    with open(path_reactions, "rb") as file_source:
        reactions = pickle.load(file_source)
    with open(path_metabolites, "rb") as file_source:
        metabolites = pickle.load(file_source)
    filtration_compartments = utility.read_file_table(
        path_file=path_filtration_compartments,
        names=None,
        delimiter="\t"
    )
    filtration_processes = utility.read_file_table(
        path_file=path_filtration_processes,
        names=None,
        delimiter="\t"
    )
    simplification_reactions = utility.read_file_table(
        path_file=path_simplification_reactions,
        names=None,
        delimiter="\t"
    )
    simplification_metabolites = utility.read_file_table(
        path_file=path_simplification_metabolites,
        names=None,
        delimiter="\t"
    )
    # Compile and return information.
    return {
        "compartments": compartments,
        "processes": processes,
        "reactions": reactions,
        "metabolites": metabolites,
        "filtration_compartments": filtration_compartments,
        "filtration_processes": filtration_processes,
        "simplification_reactions": simplification_reactions,
        "simplification_metabolites": simplification_metabolites
    }


def collect_candidate_reactions(
    reactions=None,
    filtration_compartments=None,
    filtration_processes=None,
    simplification_reactions=None,
    simplification_metabolites=None
):
    """
    Collects information about candidate reactions.

    Evaluate the candidacy of each reaction.
    A reaction's candidacy depends on the relevance of its behavior in the
    context of the relevance of compartmentalization in general and the
    relevance of specific compartments and processes.
    A reaction's candidacy also depends on the relevance of specific
    metabolites that are its participants.

    Factors in a reaction's candidacy:
    1. compartmentalization
    2. filters by compartments
    3. filters by processes
    4. simplification of metabolites
    5. simplification of reactions

    arguments:
        reactions (dict<dict>): information about reactions
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        filtration_processes (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific processes
        simplification_reactions (list<dict<str>>): information about whether
            to simplify representations of specific reactions
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (dict<dict>): information about candidate reactions

    """

    pass


###############################################################################
# Procedure


def execute_procedure(directory=None):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to evaluate the candidacy of metabolites
    and reactions for representation in a network.

    arguments:
        directory (str): path to directory for source and product files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(directory=directory)
    # TODO: 2. iterate over those candidate reactions to collect candidate metabolites
    # TODO: note that at this point it is appropriate to either omit or replicate
    # TODO: a reaction's metabolites according to the simplification method
    # TODO: ie consider simplification of metabolites when deciding candidacy
    # TODO: 3. the network procedure should not need to consider filtration
    # TODO: or simplification at all. It just creates nodes and links and transfers
    # TODO: over info from the candidates and metabolites and reactions
    # Evaluate reactions' relevance.
    # Collect candidate reactions.
    reactions_candidates = collect_candidate_reactions(
        reactions=source["reactions"],
        filtration_compartments=source["filtration_compartments"],
        filtration_processes=source["filtration_processes"],
        simplification_reactions=source["simplification_reactions"],
        simplification_metabolites=source["simplification_metabolites"]
    )
    # Collect candidate metabolites.
    metabolites_candidates = collect_candidate_metabolites(
        metabolites=source["metabolites"],
        compartments=source["compartments"],
        reactions_candidates=reactions_candidates,
        filtration_compartments=source["filtration_compartments"],
        filtration_processes=source["filtration_processes"],
        simplification_metabolites=source["simplification_metabolites"]
    )
    # TODO: Prepare some sort of report of candidate reactions and metabolites.
    # TODO: maybe do degrees (counts reactions/metabolites) to inform simplification



    pass
