"""
Author:

    Thomas Cameron Waller
    tcameronwaller@gmail.com
    Department of Biochemistry
    University of Utah
    Room 4100, Emma Eccles Jones Medical Research Building
    15 North Medical Drive East
    Salt Lake City, Utah 84112
    United States of America

License:

    This file is part of MetaboNet
    (https://github.com/tcameronwaller/metabonet/).

    MetaboNet supports definition and analysis of custom metabolic networks.
    Copyright (C) 2019 Thomas Cameron Waller

    MetaboNet is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    MetaboNet is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with MetaboNet. If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard.
import setuptools

# Relevant.

# Custom.

#dir()
#importlib.reload()

###############################################################################
# Functionality

###############################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    setuptools.setup(
        name="metabonet",
        version="0.0.1",
        description="Curation, Definition, and Analysis of Metabolic Networks",
        author="T. Cameron Waller",
        author_email="tcameronwaller@gmail.com",
        url="https://github.com/tcameronwaller/metabonet",
        packages=[
            "metabonet", "metabonet.metabocurator"
        ],#setuptools.find_packages(),
        #package_dir={"":"metabonet"},
        install_requires=[
            "pandas",
            "numpy",
            "networkx",
            "matplotlib",
            "wordcloud",
        ],
        license="https://www.gnu.org/licenses/gpl.html",
        entry_points={
            "console_scripts": [
                "metabonet = metabonet.interface:execute_procedure"
            ]
        }
    )



if (__name__ == "__main__"):
    execute_procedure()
