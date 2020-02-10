## Installation

    conda create -n charger_py38 python=3.8 poetry pip
    poetry install

## Usage

    charger -h

## Development
Run all the test by:

    pytest

The repo use isort to sort the import order, black to format the code, flake8 to check coding style, and mypy to check typing. Run them before commit:

    isort
    black src tests
    flake8
    mypy src

The repo comes with VSCode settings. In VSCode, black, mypy and flake8 will be run on every file save.

Build the documentation using Sphinx:

    cd docs
    make html

And the documentation will be available under `docs/_build/html`.