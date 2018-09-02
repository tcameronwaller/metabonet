"""
"""

###############################################################################
# Notes


###############################################################################
# Installation and importation of packages and modules

# Packages and modules from the python standard library

import os
#import sys
import shutil
#import importlib
import csv
import copy
import pickle

# Packages and modules from third parties

#import numpy
#import pandas

# Packages and modules from local source

#dir()
#importlib.reload()

###############################################################################
# Functionality


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


def copy_interpret_content_recon2m2(content=None):
    """
    Copies and interprets content from Recon 2M.2

    This function copies and interprets content from a metabolic model in
    Systems Biology Markup Language (SBML), a form of Extensible Markup
    Language (XML).

    arguments:
        content (object): content from Recon 2M.2 in SBML

    raises:

    returns:
        (object): references to definition of name space and sections within
            content

    """

    # Copy content.
    content_copy = copy.deepcopy(content)
    # Define name space.
    space = {
        "version": "http://www.sbml.org/sbml/level2/version4",
        "syntax": "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    }
    # Set references to sections within content.
    sbml = content_copy.getroot()
    model = sbml[0]
    compartments = model[1]
    metabolites = model[2]
    reactions = model[3]
    # Compile and return information.
    return {
        "space": space,
        "content": content_copy,
        "model": model,
        "compartments": compartments,
        "metabolites": metabolites,
        "reactions": reactions
    }


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
    with open(path_file, "w", newline="") as file_product:
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

    returns:
        (list): values from records

    raises:

    """

    def access(record):
        return record[key]
    return list(map(access, records))


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


def interpret_content_hmdb(content=None):
    """
    Copies and interprets content from Human Metabolome Database (HMDB)

    This function copies and interprets content from the Human Metabolome
    Database in Extensible Markup Language (XML).

    arguments:
        content (object): content from HMDB in XML

    returns:
        (object): references to definition of name space and sections within
            content

    raises:

    """

    # Define name space
    space = {
        "base": "http://www.hmdb.ca"
    }
    # Set references to sections within content
    metabolites = content.getroot()
    #for child in metabolites[0]:
    #    print(child.tag)
    # Compile information
    return {
        "space": space,
        "content": content,
        "metabolites": metabolites
    }


###############################################################################
# Procedure
