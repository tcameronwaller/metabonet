"""
Module to provide template of common structure of modules.

Title:

    utility

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
import csv
import copy

# Relevant

# Custom

#dir()
#importlib.reload()

###############################################################################
# Functionality


# General.


def confirm_path_directory(path=None):
    """
    Confirms that a path to a directory exists.

    Creates a directory if it does not already exist.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if not os.path.exists(path):
        os.makedirs(path)


def read_file_table(path_file=None, names=None, delimiter=None):
    """
    Reads and organizes source information from file

    This function reads and organizes relevant information from file.

    arguments:
        path_file (str): path to directory and file
        names (list<str>): names for values in each row of table
        delimiter (str): delimiter between values in the table

    returns:
        (list<dict>): tabular information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        reader = csv.DictReader(
            file_source, fieldnames=names, delimiter=delimiter
        )
        information = list(map(lambda row: dict(row), list(reader)))
    # Return information
    return information


def write_file_table(
    information=None, path_file=None, names=None, delimiter=None
):
    """
    Writes information to file

    arguments:
        path_file (str): path to directory and file
        names (list<str>): names for values in each row of table
        delimiter (str): delimiter between values in the table

    returns:

    raises:

    """

    # Write information to file
    #with open(out_file_path_model, "w") as out_file:
    #    out_file.write(content_identifier)
    with open(path_file, "w") as file_product:
        writer = csv.DictWriter(
            file_product, fieldnames=names, delimiter=delimiter
        )
        writer.writeheader()
        writer.writerows(information)


def find(match=None, sequence=None):
    """
    Finds the first element in a sequence to match a condition, otherwise none

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (object | NoneType): first element from sequence to match condition or
            none

    raises:

    """

    for element in sequence:
        if match(element):
            return element
    return None


def find_index(match=None, sequence=None):
    """
    Finds index of first element in sequence to match a condition, otherwise -1

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (int): index of element if it exists or -1 if it does not exist

    raises:

    """

    for index, element in enumerate(sequence):
        if match(element):
            # Element matches condition
            # Return element's index
            return index
    # Not any elements match condition
    # Return -1
    return -1


def find_all(match=None, sequence=None):
    """
    Finds all elements in a sequence to match a condition, otherwise none.

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (list<dict> | NoneType): elements from sequence to match condition or
            none

    raises:

    """

    matches = []
    for element in sequence:
        if match(element):
            matches.append(element)
    if len(matches) > 0:
        return matches
    else:
        return None


def collect_unique_elements(elements_original=None):
    """
    Collects unique elements

    arguments:
        elements_original (list): sequence of elements

    returns:
        (list): unique elements

    raises:

    """

    elements_novel = []
    for element in elements_original:
        if element not in elements_novel:
            elements_novel.append(element)
    return elements_novel


def collect_value_from_records(key=None, records=None):
    """
    Collects a single value from multiple records

    arguments:
        key (str): key of value in each record
        records (list<dict>): sequence of records

    raises:

    returns:
        (list): values from records

    """

    def access(record):
        return record[key]
    return list(map(access, records))


def collect_values_from_records(key=None, records=None):
    """
    Collects values from multiple records.

    arguments:
        key (str): key of value in each record
        records (list<dict>): sequence of records

    raises:

    returns:
        (list): values from records

    """

    collection = []
    for record in records:
        collection.extend(record[key])
    return collection


def compare_lists_by_inclusion(list_one=None, list_two=None):
    """
    Compares lists by inclusion

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (bool): whether first list includes all elements from second

    raises:

    """

    def match(element_two=None):
        return element_two in list_one
    matches = list(map(match, list_two))
    return all(matches)


def compare_lists_by_mutual_inclusion(list_one=None, list_two=None):
    """
    Compares lists by mutual inclusion

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (bool): whether each list includes all elements from the other

    raises:

    """

    forward = compare_lists_by_inclusion(
        list_one=list_one,
        list_two=list_two
    )
    reverse = compare_lists_by_inclusion(
        list_one=list_two,
        list_two=list_one
    )
    return forward and reverse


def filter_common_elements(list_one=None, list_two=None):
    """
    Filters elements by whether both of two lists include them

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (list): elements that both of two lists include

    raises:

    """

    def match(element_two=None):
        return element_two in list_one
    return list(filter(match, list_two))


def collect_records_targets_by_categories(
    target=None,
    category=None,
    records=None
):
    """
    Collects values of a target attribute that occur together in records with
    each value of another category attribute.
    Each record has a single value of the target attribute.
    Each record can have either a single value or multiple values of the
    category attribute.
    These collections do not necessarily include only unique values of the
    target attribute.

    arguments:
        target (str): name of attribute in records to collect for each category
        category (str): name of attribute in records to define categories
        records (list<dict>): records with target and category attributes

    raises:

    returns:
        (dict<list<str>>): values of the target attribute that occur together
            in records with each value of the category attribute

    """

    def collect_record_target_by_category(
        target_value=None,
        category_value=None,
        collection_original=None,
    ):
        collection_novel = copy.deepcopy(collection_original)
        # Determine whether collection includes the category's value.
        if category_value in collection_novel.keys():
            # Collection includes the category's value.
            target_values = collection_novel[category_value]
            target_values.append(target_value)
            # Include target's value in collection.
            collection_novel[category_value] = target_values
        else:
            # Collection does not include the category's value.
            # Include category's value and target's value in collection.
            collection_novel[category_value] = [target_value]
        return collection_novel

    collection = {}
    for record in records:
        target_value = record[target]
        category_values = record[category]
        if isinstance(category_values, list):
            for category_value in category_values:
                collection = collect_record_target_by_category(
                    target_value=target_value,
                    category_value=category_value,
                    collection_original=collection
                )
        else:
            category_value = category_values
            collection = collect_record_target_by_category(
                target_value=target_value,
                category_value=category_value,
                collection_original=collection
            )
    return collection


def collect_values_from_records_in_reference(
    key=None, identifiers=None, reference=None
):
    """
    Collects a single value from a specific record in a reference.

    arguments:
        key (str): key of value in record
        identifiers (list<str>): identifiers of records in reference
        reference (dict<dict<str>>): reference of records

    raises:

    returns:
        (list<str>): values from records

    """

    values = []
    for identifier in identifiers:
        record = reference[identifier]
        value = record[key]
        values.append(value)
    return values


def filter_nonempty_elements(elements_original=None):
    """
    Filters nonempty elements.

    arguments:
        elements_original (list<str>): sequence of elements

    returns:
        (list<str>): non-empty elements

    raises:

    """

    elements_novel = []
    for element in elements_original:
        if len(element) > 0:
            elements_novel.append(element)
    return elements_novel


# Human Metabolome Database (HMDB).


def match_hmdb_entries_by_identifiers_names(
    identifiers=None,
    names=None,
    summary_hmdb=None
):

    """
    Matches entries from Human Metabolome Database by identifiers or names.

    arguments:
        identifiers (list<str>): identifiers by which to find entries in HMDB
        names (list<str>): names by which to find entries in HMDB
        summary_hmdb (dict<dict>): information about metabolites from Human
            Metabolome Database (HMDB)

    returns:
        (list<str>): keys of entries in HMDB

    raises:

    """

    # Test.
    #hmdb_keys = utility.match_hmdb_entries_by_identifiers_names(
    #    identifiers=[],
    #    names=["pyruvate"],
    #    summary_hmdb=source["summary_hmdb"]
    #)
    # HMDB0000243

    # Ensure identifiers and names are not empty.
    identifiers_valid = filter_nonempty_elements(identifiers)
    names_valid = filter_nonempty_elements(names)
    # Determine whether measurement's record include reference to HMDB.
    if (len(identifiers_valid) > 0):
        # Measurement's record includes references to HMDB.
        # Match measurement's record to a entries in HMDB.
        # Match by identifier.
        hmdb_keys = filter_hmdb_entries_by_identifiers(
            identifiers=identifiers_valid,
            summary_hmdb=summary_hmdb
        )
    elif (len(names_valid) > 0):
        # Measurement's record does not include reference to HMDB.
        # Match measurement's record to an entry in HMDB.
        # Attempt to match by name.
        hmdb_keys = filter_hmdb_entries_by_synonyms(
            names=names_valid,
            summary_hmdb=summary_hmdb
        )
    else:
        hmdb_keys = []
    # Return information.
    return hmdb_keys


def filter_hmdb_entries_by_identifiers(
    identifiers=None, summary_hmdb=None
):
    """
    Filters entries from HMDB by their identifiers.

    arguments:
        identifiers (list<str>): identifiers by which to find entries in HMDB
        summary_hmdb (dict<dict>): information about metabolites from Human
            Metabolome Database (HMDB)

    returns:
        (list<str>): keys of entries in HMDB

    raises:

    """

    keys = []
    for key, record in summary_hmdb.items():
        hmdb_entry_identifiers = record["references_hmdb"]
        # Determine whether any of entry's identifiers match the metabolite's
        # references
        checks = []
        for identifier in identifiers:
            check = identifier in hmdb_entry_identifiers
            checks.append(check)
        if any(checks):
            # The entry matches the metabolite's references.
            keys.append(key)
    return keys


def filter_hmdb_entries_by_synonyms(
    names=None, summary_hmdb=None
):
    """
    Filters entries from HMDB by their synonyms.

    arguments:
        names (list<str>): names by which to find entries in HMDB
        summary_hmdb (dict<dict>): information about metabolites from Human
            Metabolome Database (HMDB)

    returns:
        (list<str>): keys of entries in HMDB

    raises:

    """

    keys = []
    for key, record in summary_hmdb.items():
        synonyms = record["synonyms"]
        synonyms_lower = []
        for synonym in synonyms:
            synonyms_lower.append(synonym.lower())
        # Determine whether any of entry's identifiers match the metabolite's
        # references
        checks = []
        for name in names:
            name_lower = name.lower()
            check = name_lower in synonyms_lower
            checks.append(check)
        if any(checks):
            # The entry matches the metabolite's references
            keys.append(key)
    return keys


# Metabolic information.


def collect_reaction_participants_value(
    key=None, criteria=None, participants=None
):
    """
    Collects a value from a reaction's specific participants

    arguments:
        key (str): key of value to collect from each participant
        criteria (dict<list>): criteria by which to select participants
        participants (list<dict>): information about a reaction's participants

    returns:
        (list<str>): values from a reaction's participants

    raises:

    """

    participants_match = filter_reaction_participants(
        criteria=criteria, participants=participants
    )
    return collect_value_from_records(key=key, records=participants_match)


def filter_reaction_participants(criteria=None, participants=None):
    """
    Filters a reaction's participants by multiple criteria

    arguments:
        criteria (dict<list>): criteria by which to select participants
        participants (list<dict>): information about a reaction's participants

    returns:
        (list<dict>): information about a reaction's participants

    raises:

    """

    def match(participant):
        if "metabolites" in criteria:
            match_metabolite = (
                participant["metabolite"] in criteria["metabolites"]
            )
        else:
            match_metabolite = True
        if "compartments" in criteria:
            match_compartment = (
                participant["compartment"] in criteria["compartments"]
            )
        else:
            match_compartment = True
        if "roles" in criteria:
            match_role = participant["role"] in criteria["roles"]
        else:
            match_role = True
        return match_metabolite and match_compartment and match_role
    return list(filter(match, participants))


###############################################################################
# Procedure
