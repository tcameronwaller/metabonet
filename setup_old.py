"""
Module to accept parameters from user and execute procedures.

Title:

    metabonet

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
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

# Standard.
import os
import sys
import subprocess
import setuptools
from setuptools.command.develop import develop
from setuptools.command.install import install

# Relevant.

# Custom.
def test_helps():

    return_list = []

    # Check help menus
    return_list.append(
        subprocess.call(
            'metabonet --help > /dev/null',
            shell = True
            )
        )
    return_list.append(
        subprocess.call(
            'metabonet model --help > /dev/null',
            shell = True
            )
        )
    return_list.append(
        subprocess.call(
            'metabonet network --help > /dev/null',
            shell = True
            )
        )
    return_list.append(
        subprocess.call(
            'metabonet clean --help > /dev/null',
            shell = True
            )
        )

    # Ensure all outputs were successful
    result = all(x == 0 for x in return_list)
    if result == False:
        subprocess.call(
            'pip uninstall metabonet -y',
            shell = True
            )
        subprocess.call(
            'echo "Error: Something went wrong during installation. Please fix the errors and try reinstalling"',
            shell = True
            )
        sys.exit(1)
    else:
        subprocess.call(
            'echo "Tests passed without error, installation successful"',
            shell = True
            )

# Install Recon and HMDB databases to dock/source/ location in metabonet
def install_databases():

    __path__  =  os.path.dirname(os.path.realpath(__file__))
    __path__ = __path__ + '/dock/source/'

    subprocess.call(
        ('echo "Installing databases..."; '
        + 'curl https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml?download=1 -o ' + str(__path__) + 'recon2m2.xml; '
        + 'echo "Recon2M.2_MNX_Entrex_Gene Human Metabolism database installed"; '
        + 'curl http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip -o ' + str(__path__) + 'hmdb_metabolites.zip; '
        + 'unzip ' + str(__path__) + 'hmdb_metabolites.zip -d ' + str(__path__) + '; '
        + 'rm ' + str(__path__) + 'hmdb_metabolites.zip; '
        + 'echo "HMDB metabolite database v4.0 installed"; '
        + 'echo "Database installs complete"'),
        shell = True
        )

class PostDevelopCommand(develop):
    # Post-installation for python setup.py develop
    def run(self):
        develop.run(self)
        test_helps()
        install_databases()

class PostInstallCommand(install):
    # Post-installation for python setup.py install
    def run(self):
        install.run(self)
        test_helps()
        install_databases()

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
              'pandas',
              'numpy',
              'networkx',
              'wordcloud'
              ],
        license="https://www.gnu.org/licenses/gpl.html",
        cmdclass={
            'develop': PostDevelopCommand,
            'install': PostInstallCommand
            },
        entry_points={
            "console_scripts": [
                "metabonet = metabonet.interface:execute_procedure"
            ]
        }
    )


if (__name__ == "__main__"):
    execute_procedure()
