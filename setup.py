#!/usr/bin/env python
# -*- coding: utf-8 -*-

# NOTE: code adapted from https://github.com/navdeep-G/setup.py/blob/master/setup.py

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

from configparser import ConfigParser
import io
import os
from pathlib import Path
import sys
from shutil import rmtree
from setuptools import find_packages, setup, Command
from configparser import ConfigParser
from pathlib import Path

def clean(s):
    return s.strip().strip("'\"")

def clean(s):
    return s.strip().strip("'\"")

# Package meta-data.
<<<<<<< HEAD
=======
# NAME = 'ncbi-submit'
# DESCRIPTION = 'A tool for submitting to NCBI (SRA, BioSample, & GenBank).'
# URL = 'https://github.com/enviro-lab/ncbi_interact'
# EMAIL = 'skunklem@uncc.edu'
# AUTHOR = 'Sam Kunkleman'
# REQUIRES_PYTHON = '>=3.8'
# VERSION = ''
>>>>>>> ed450ef62c4e18b969e601f60621894e67ab6b7e
toml = Path(__file__).resolve().parent/"pyproject.toml"
config = ConfigParser(converters={
    'list': lambda x: [clean(i) for i in x.strip('[]').split(',')],
    'clean': lambda x: clean(x)
})
config.read(toml)
authors = config.getlist("tool.poetry","authors")
NAME = config.getclean("tool.poetry","name")
DESCRIPTION = config.getclean("tool.poetry","description")
URL = config.getclean("tool.poetry","repository")
AUTHOR,EMAIL = authors[0].strip(">").split(" <")
REQUIRES_PYTHON = config.getclean("tool.poetry.dependencies","python")
VERSION = config.getclean("tool.poetry","version")
LICENSE = config.getclean("tool.poetry","license")

# What packages are required for this module to be executed?
<<<<<<< HEAD
=======
# REQUIRED = [
#     'biopython',
#     'pandas',
# ]
>>>>>>> ed450ef62c4e18b969e601f60621894e67ab6b7e
REQUIRED = [f"{k}{clean(v)}" for k,v in config["tool.poetry.dependencies"].items() if "python" not in k]

# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, 'version.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
<<<<<<< HEAD
    packages=find_packages(
        # exclude=["tests", "*.tests", "*.tests.*", "tests.*"]
    ),
=======
    # packages=find_packages(exclude=["ncbi_submit", "example"]),
    packages=find_packages(
        # exclude=["tests", "*.tests", "*.tests.*", "tests.*"]
    ),
    # packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
>>>>>>> ed450ef62c4e18b969e601f60621894e67ab6b7e
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['ncbi_submit'],

    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ],
    # # $ setup.py publish support.
    # cmdclass={
    #     'upload': UploadCommand,
    # },
)