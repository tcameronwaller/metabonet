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
    compartmentalization=None,
    filtration_compartments=None,
    filtration_processes=None,
    simplification_reactions=None,
    simplification_metabolites=None
):
    """
    Collects information about candidate reactions.

    arguments:
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
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

    reactions_candidacy = {}
    for reaction in reactions.values():
        match, record = collect_candidate_reaction(
            reaction_identifier= reaction["identifier"],
            reactions=reactions,
            compartmentalization=compartmentalization,
            filtration_compartments=filtration_compartments,
            filtration_processes=filtration_processes,
            simplification_reactions=simplification_reactions,
            simplification_metabolites=simplification_metabolites
        )
        if match:
            reactions_candidacy[record["identifier"]] = record
    return reactions_candidacy


def collect_candidate_reaction(
    reaction_identifier=None,
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    filtration_processes=None,
    simplification_reactions=None,
    simplification_metabolites=None
):
    """
    Collects information about a candidate reaction.

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
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
        (tuple<bool, dict>): whether reaction is a candidate, and information
            about the reaction candidate

    """

    # Determine reaction's candidacy.
    candidacy = determine_reaction_candidacy(
        reaction_identifier= reaction_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        filtration_processes=filtration_processes,
        simplification_reactions=simplification_reactions,
        simplification_metabolites=simplification_metabolites
    )


def determine_reaction_candidacy(
    reaction_identifier=None,
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    filtration_processes=None,
    simplification_reactions=None,
    simplification_metabolites=None
):
    """
    Determines whether a reaction is a candidate for representation in a
    network.

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
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
        (bool): whether reaction is a candidate

    """

    # Relevance.
    relevance = determine_reaction_relevance(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        filtration_processes=filtration_processes,
        simplification_reactions=simplification_reactions,
        simplification_metabolites=simplification_metabolites
    )
    # Redundancy.
    redundancy = determine_reaction_redundancy(
        reaction_identifier= reaction_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        filtration_processes=filtration_processes,
        simplification_reactions=simplification_reactions,
        simplification_metabolites=simplification_metabolites
    )
    # Determine whether reaction is a candidate.
    return relevance and not redundancy


def determine_reaction_relevance(
    reaction_identifier=None,
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    filtration_processes=None,
    simplification_reactions=None,
    simplification_metabolites=None
):
    """
    Determines whether a reaction is relevant.

    A reaction's relevance depends on the relevance of its behavior in the
    context of the relevance of compartmentalization in general and the
    relevance of specific compartments and processes.
    A reaction's relevance also depends on the relevance of specific
    metabolites that are its participants.

    Factors in a reaction's relevance:
    1. compartmentalization
    2. filters by compartments
    3. filters by processes
    4. simplification of metabolites
    5. simplification of reactions

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
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
        (bool): whether reaction is relevant

    """

    # Simplification.
    simplification = determine_reaction_simplification(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        simplification_reactions=simplification_reactions
    )
    # Process.
    process = determine_reaction_process_relevance(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        filtration_processes=filtration_processes
    )
    # Behavior.
    behavior = determine_reaction_behavior_relevance(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        simplification_metabolites=simplification_metabolites
    )
    # Determine whether reaction is relevant.
    return (not simplification) and process and behavior


def determine_reaction_simplification(
    reaction_identifier=None,
    simplification_reactions=None
):
    """
    Determines whether a reaction has designation for simplification.

    arguments:
        reaction_identifier (str): identifier of a reaction
        simplification_reactions (list<dict<str>>): information about whether
            to simplify representations of specific reactions

    raises:

    returns:
        (bool): whether reaction has designation for simplification

    """

    def match_reaction_simplification(record):
        return record["identifier"] == reaction_identifier
    match = utility.find(
        match=match_reaction_simplification,
        sequence=simplification_reactions
    )
    return match is not None


def determine_reaction_process_relevance(
    reaction_identifier=None,
    reactions=None,
    filtration_processes=None
):
    """
    Determines whether a reaction belongs to any relevant processes.

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        filtration_processes (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific processes

    raises:

    returns:
        (bool): whether reaction belongs to relevant process

    """

    reaction = reactions[reaction_identifier]
    reaction_processes = reaction["processes"]
    qualifiers = []
    for process in reaction_processes:
        relevance = determine_set_relevance(
            identifier=process,
            filters=filtration_processes
        )
        if relevance:
            qualifiers.append(process)
    return len(qualifiers) > 0


def determine_set_relevance(identifier=None, filters=None):
    """
    Determines whether a metabolic set, compartment or process, is relevant.

    arguments:
        identifier (str): identifier of a set
        filters (list<dict<str>>): information about relevance of specific sets

    raises:

    returns:
        (bool): whether set is relevant

    """

    def match(record):
        return record["identifier"] == identifier
    record = utility.find(
        match=match,
        sequence=filters
    )
    if record is not None:
        return record["relevance"]
    else:
        # If a set does not have a record in filters, then assume that it is
        # relevant.
        return True


def determine_reaction_behavior_relevance(
    reaction_identifier=None,
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    simplification_metabolites=None
):
    """
    Determines whether a reaction's behavior is relevant.

    The relevance of a reaction's behavior depends on the relevance of
    compartmentalization in general, and the relevance of specific compartments
    and metabolites.

    Factors in relevance of a reaction's behavior:
    1. reaction's primary behavior, whether conversion or transport
    2. compartmentalization
    3. filters by compartments
    4. simplification of metabolites

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (bool): whether reaction is a candidate

    """

    reaction = reactions[reaction_identifier]
    # Determine reaction's behavior.
    if reaction.conversion:
        # Reaction involves chemical conversion.
        # Reaction's relevance requires a relevant reactant participant and
        # product participant that are chemically distinct.
        participation = determine_reaction_conversion_participation(
            reaction_identifier=reaction_identifier,
            reactions=reactions,
            filtration_compartments=filtration_compartments,
            simplification_metabolites=simplification_metabolites
        )
    elif reaction.transport:
        # Reaction does not involve chemical conversion.
        # Reaction involves compartmental transport.
        # Reaction's relevance requires relevance of compartmentalization.
        if compartmentalization:
            # Reaction's relevance requires a relevant reactant participant and
            # product participant that are chemically identical and in separate
            # compartments.
            participation = determine_reaction_transport_participation(
                reaction_identifier=reaction_identifier,
                reactions=reactions,
                filtration_compartments=filtration_compartments,
                simplification_metabolites=simplification_metabolites
            )
        else:
            participation = False
    return participation


def determine_reaction_conversion_participation(
    reaction_identifier=None,
    reactions=None,
    filtration_compartments=None,
    simplification_metabolites=None
):
    """
    Determines whether a reaction involves relevant participation for
    conversion behavior.

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (bool): whether reaction involves relevant participation

    """

    # Determine relevant participants.
    participants = determine_reaction_relevant_participants(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification_metabolites=simplification_metabolites
    )
    # Determine whether any reactant participants and product participants are
    # relevant.
    metabolites_reactant = collect_reaction_participants_value(
        key="metabolite",
        criteria={"roles": ["reactant"]},
        participants=participants
    )
    metabolites_product = collect_reaction_participants_value(
        key="metabolite",
        criteria={"roles": ["product"]},
        participants=participants
    )
    return (len(metabolites_reactant) > 0) and (len(metabolites_product > 0)


def determine_reaction_transport_participation(
    reaction_identifier=None,
    reactions=None,
    filtration_compartments=None,
    simplification_metabolites=None
):
    """
    Determines whether a reaction involves relevant participation for transport
    behavior.

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (bool): whether reaction involves relevant participation

    """

    # Determine relevant participants.
    participants = determine_reaction_relevant_participants(
        reaction_identifier=reaction_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification_metabolites=simplification_metabolites
    )
    # Determine whether any reactant participants and product participants
    # match transport.
    reaction = reactions[reaction_identifier]
    transports_original = reaction["transports"]
    transports_novel = []
    for transport in transports_original:
        metabolite = transport["metabolite"]
        compartments = transport["compartments"]
        compartments_reactant = collect_reaction_participants_value(
            key="compartment",
            criteria={
                "metabolites": [metabolite],
                "compartments": compartments,
                "roles": ["reactant"]
            },
            participants=participants
        )
        compartments_product = collect_reaction_participants_value(
            key="compartment",
            criteria={
                "metabolites": [metabolite],
                "compartments": compartments,
                "roles": ["product"]
            },
            participants=participants
        )
        if (len(compartments_reactant) > 0 and len(compartments_product) > 0):
            compartments_difference = not utility
            .compare_lists_by_mutual_inclusion(
                list_one=compartments_reactant,
                list_two=compartments_product
            )
            if compartments_difference:
                transports_novel.append(transport)
    return len(transports_novel) > 0


# TODO: I'll need a separate function for candidate participants...
# TODO: candidate participants are different from relevant participants
# TODO: candidate participants include simplification metabolites that are replicated for reaction
# TODO: candidate participants do not include omitted metabolites

def determine_reaction_relevant_participants(
    reaction_identifier=None,
    reactions=None,
    filtration_compartments=None,
    simplification_metabolites=None
):
    """
    Determines a reaction's relevant participants.

    A participant's relevance depends on the relevance of its metabolite and
    compartment.

    Factors in relevance of a reaction's participant:
    1. simplification of metabolites
    2. filters by compartments

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (list<dict<str>>): reaction's relevant participants

    """

    reaction = reactions[reaction_identifier]
    participants_original = reaction["participants"]
    # Collect relevant participants.
    participants_novel = []
    for participant in participants_original:
        # Determine relevance of participant's compartment.
        compartment = participant["compartment"]
        compartment_relevance = determine_set_relevance(
            identifier=compartment,
            filters=filtration_compartments
        )
        # Determine relevance of participant's metabolite.
        metabolite = participant["metabolite"]
        simplification = determine_metabolite_simplification(
            metabolite_identifier=metabolite,
            compartment_identifier=compartment,
            simplification_metabolites=simplification_metabolites
        )
        metabolite_relevance = (
            not simplification["omission"] and
            not simplification["replication"]
        )
        # Determine whether participant is relevant.
        if metabolite_relevance and compartment_relevance:
            participants_novel.append(participant)
    return participants_novel


def determine_metabolite_simplification(
    metabolite_identifier=None,
    compartment_identifier=None,
    simplification_metabolites=None
):
    """
    Determines a metabolite's a designation for simplification.

    A simplification must match both the metabolite and the compartment.

    arguments:
        metabolite_identifier (str): identifier of a metabolite
        compartment_identifier (str): identifier of a compartment
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (dict<bool>): whether metabolite has designation for simplification by
            omission or replication

    """

    def match_metabolite(record):
        return record["metabolite"] == metabolite_identifier
    simplifications_metabolite = utility.find_all(
        match=match_metabolite,
        sequence=simplification_metabolites
    )
    if simplifications_metabolite is not None:
        # Determine whether any of metabolite's simplifications match the
        # compartment.
        def match_compartment(record):
            return (
                (record["compartment"] == "all") or
                (record["compartment"] == compartment_identifier)
            )
        simplification_compartment = utility.find(
            match=match_compartment,
            sequence=simplifications_metabolite
        )
        if simplification_compartment is not None:
            record = {
                "omission": simplification_compartment["omission"],
                "replication": simplification_compartment["replication"]
            }
        else:
            record = {
                "omission": False,
                "replication": False
            }
    else:
        record = {
            "omission": False,
            "replication": False
        }
    return record


def determine_reaction_redundancy(
    reaction_identifier=None,
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    filtration_processes=None,
    simplification_reactions=None,
    simplification_metabolites=None
):
    """
    Determines whether a reaction is redundant.

    A reaction is redundant if its participants are identical to those of
    another relevant reaction.
    The relevance of compartmentalization determines whether compartments
    influence this comparison of participants.

    arguments:
        reaction_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
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
        (bool): whether reaction is redundant

    """

    reaction = reactions["reaction_identifier"]
    replicates = reaction["replicates"]
    # Determine relevant replicates.
    def match_relevance(replicate):
        identity = (replicate == reaction_identifier)
        relevance = determine_reaction_relevance(
            reaction_identifier=replicate,
            reactions=reactions,
            compartmentalization=compartmentalization,
            filtration_compartments=filtration_compartments,
            filtration_processes=filtration_processes,
            simplification_reactions=simplification_reactions,
            simplification_metabolites=simplification_metabolites
        )
        return not identity and relevance
    replicates_relevant = list(filter(match_relevance, replicates))
    # Determine redundant replicates.
    def match_redundancy(replicate):
        redundancy = determine_replicate_reactions_redundancy(
            reaction_one_identifier=reaction_identifier,
            reaction_two_identifier=replicate,
            reactions=reactions,
            compartmentalization=compartmentalization,
            filtration_compartments=filtration_compartments,
            simplification_metabolites=simplification_metabolites
        )
        return not redundancy
    replicates_redundant = list(filter(match_redundancy, replicates_relevant))
    # Determine whether reaction is priority replicate.
    if len(replicates_redundant) > 0:
        # TODO: next figure out priority reaction from relevant, redundant replicates
        # TODO: need to rank the reactions...



    # Determine whether replicate reaction is redundant.
    # Consider reversibility and compartmentalization.

    # If any replicates are both relevant and redundant, then determine whether
    # reaction is the proirity.
    # Use a simple ranking algorithm to compare reactions.



def determine_replicate_reactions_redundancy(
    reaction_one_identifier=None,
    reaction_two_identifier=None
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    simplification_metabolites=None
):
    """
    Determines whether two replicate reactions are redundant.

    arguments:
        reaction_one_identifier (str): identifier of a reaction
        reaction_two_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (bool): whether reactions are redundant

    """

    reaction_one = reactions[reaction_one_identifier]
    reaction_two = reactions[reaction_two_identifier]
    # Determine whether reactions have redundant reversibility.
    reversibility = (
        reaction_one["reversibility"] == reaction_two["reversibility"]
    )
    # Determine whether reactions have redundant participants.
    # Determine relevant participants.
    participation = determine_reactions_participation_redundancy(
        reaction_one_identifier=reaction_one_identifier,
        reaction_two_identifier=reaction_two_identifier,
        reactions=reactions,
        compartmentalization=compartmentalization,
        filtration_compartments=filtration_compartments,
        simplification_metabolites=simplification_metabolites
    )
    # Determine whether reactions are redundant.
    return reversibility and participation


def determine_reactions_participation_redundancy(
    reaction_one_identifier=None,
    reaction_two_identifier=None
    reactions=None,
    compartmentalization=None,
    filtration_compartments=None,
    simplification_metabolites=None
):
    """
    Determines whether two replicate reactions have redundant participation.

    arguments:
        reaction_one_identifier (str): identifier of a reaction
        reaction_two_identifier (str): identifier of a reaction
        reactions (dict<dict>): information about reactions
        compartmentalization (bool): whether compartmentalization is relevant
        filtration_compartments (list<dict<str>>): information about whether to
            remove metabolites and reactions relevant to specific compartments
        simplification_metabolites (list<dict<str>>): information about whether
            to simplify representations of specific metabolites

    raises:

    returns:
        (bool): whether reactions' participants are redundant

    """

    # Determine relevant participants.
    participants_one = determine_reaction_relevant_participants(
        reaction_identifier=reaction_one_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification_metabolites=simplification_metabolites
    )
    participants_two = determine_reaction_relevant_participants(
        reaction_identifier=reaction_two_identifier,
        reactions=reactions,
        filtration_compartments=filtration_compartments,
        simplification_metabolites=simplification_metabolites
    )
    # Determine whether participants are redundant.
    if compartmentalization:
        redundancy = determine_participants_attributes_mutual_redundancy(
            participants_one=participants_one,
            participants_two=participants_two,
            attributes=["metabolite", "role", "compartment"]
        )
    else:
        redundancy = determine_participants_attributes_mutual_redundancy(
            participants_one=participants_one,
            participants_two=participants_two,
            attributes=["metabolite", "role"]
        )
    return redundancy


def determine_participants_attributes_mutual_redundancy(
    participants_one=None,
    participants_two=None,
    attributes=None
):
    """
    Determines whether participants of two reactions are redundant.

    arguments:
        participants_one (str): participants of a reaction
        participants_two (str): participants of a reaction
        attributes (list<str>): names of attributes by which to compare
            participants

    raises:

    returns:
        (bool): whether reactions' participants are redundant

    """

    comparison_one = determine_participants_attributes_redundancy(
        participants_one=participants_one,
        participants_two=participants_two,
        attributes=attributes
    )
    comparison_two = determine_participants_attributes_redundancy(
        participants_one=participants_two,
        participants_two=participants_one,
        attributes=attributes
    )
    return comparison_one and comparison_two


def determine_participants_attributes_redundancy(
    participants_one=None,
    participants_two=None,
    attributes=None
):
    """
    Determines whether participants of two reactions are redundant.

    arguments:
        participants_one (str): participants of a reaction
        participants_two (str): participants of a reaction
        attributes (list<str>): names of attributes by which to compare
            participants

    raises:

    returns:
        (bool): whether reactions' participants are redundant

    """

    iteration_one = []
    for participant_one in participants_one:
        iteration_two = []
        for participant_two in participants_two:
            iteration_three = []
            for attribute in attributes:
                match = (
                    participant_one[attribute] == participant_two[attribute]
                )
                iteration_three.append(match)
            iteration_two.append(all(iteration_three))
        iteration_one.append(any(iteration_two))
    return all(iteration_one)




###############################################################################
# Procedure


def execute_procedure(
    compartmentalization=None,
    directory=None
):
    """
    Function to execute module's main behavior.

    The purpose of this procedure is to evaluate the candidacy of metabolites
    and reactions for representation in a network.

    arguments:
        compartmentalization (bool): whether compartmentalization is relevant
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
    reactions_candidacy = collect_candidate_reactions(
        reactions=source["reactions"],
        compartmentalization=compartmentalization,
        filtration_compartments=source["filtration_compartments"],
        filtration_processes=source["filtration_processes"],
        simplification_reactions=source["simplification_reactions"],
        simplification_metabolites=source["simplification_metabolites"]
    )
    # Collect candidate metabolites.
    metabolites_candidacy = collect_candidate_metabolites(
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
