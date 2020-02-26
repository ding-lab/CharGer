Development Guide
=================

Here are details of contributing to CharGer.


Workflow
--------
A typical workflow of updating CharGer involves these steps:

1. Add new code
2. Add new tests for the new code
3. Pass all the style checks (isort_, black_, and flake8_)
4. Pass all the type checks (mypy_)
5. Pass all the functional tests (pytest_)
6. Update documentation
7. Commit the change

Below are the details of each step.

.. _isort: https://github.com/timothycrosley/isort/
.. _black: https://github.com/psf/black
.. _flake8: https://flake8.pycqa.org/
.. _mypy: http://www.mypy-lang.org/
.. _pytest: https://docs.pytest.org/


Environment setup
-----------------
Make sure the conda is available and has set up bioconda_ channels.
Create a new conda environment (for example, ``charger_py37``)::

    conda create -n charger_py37 python=3.7 pip poetry cyvcf2 pysam
    conda activate charger_py37
    git clone https://github.com/ding-lab/CharGer charger
    cd charger
    poetry install      # install charger and all its dependencies
    pre-commit install  # enforce style check at every commit

On Python 3.8, exclude cyvcf2 and pysam as bioconda hasn't started to built Python 3.8 packages::

    conda create -n charger_py38 python=3.8 poetry pip numpy
    # ... continue the same as above

.. _bioconda: https://bioconda.github.io/


Style checks
------------
Style checks are enforced to pass before any git commit (by running ``pre-commit install`` once).
To run all the style checks at any time, run the following command::

    pre-commit run -a

Otherwise, run the style checks manually by::

    isort               # Sort the import order
    black src tests     # Format the code
    flake8              # Check coding style


Type checks
-----------
Type checks use mypy_ to infer the data types of the Python variables and report any potential mismatches (say, passing string where an integer is required). Run type checks by::

    mypy --pretty src tests


Functional tests
----------------
Functional tests check if the program works as expected. All the tests are under ``tests`` folder. Some tests check one particular functions, and some other tests check if the output is expected given a certain CharGer config.

Run functional tests by::

    pytest -v

Functional tests and type checks are not run automatically for every commit, but it's recommended to run and pass all the tests and type checks.

The repo should always pass all the tests (style checks, type checks, and functional tests) described above.


Test multiple Python versions
-----------------------------
To make sure CharGer works on all Python versions, we use tox to run tests on each supported Python version::

    conda install tox tox-conda
    tox -p auto -c tox_conda.ini


Build documentation
-------------------
CharGer's documentation is powered by sphinx_ under ``docs``. Build or update the documentation by::

    cd docs
    make html

And the documentation will be available under ``docs/_build/html``.

.. _sphinx: https://www.sphinx-doc.org/


Develop with Visual Studio Code
-------------------------------
Here are some additional setup that utilize `Visual Studio Code`_'s IDE:

- Run style checks (black_, mypy_ and flake8_) at every file save
- The default build task will call sphinx_ to build the documentation
- A debug shortcut to go into CharGer internals using a startup script ``scripts/debug_example.py``

Add the following workspace settings ``.vscode/settings.json``:

.. code-block:: json

    {
        "python.autoComplete.extraPaths": ["src"],
        "python.formatting.provider": "black",
        "editor.formatOnSave": true,
        "editor.wordWrapColumn": 120,
        "python.linting.enabled": true,
        "python.linting.flake8Enabled": true,
        "python.linting.mypyEnabled": true,
        "python.linting.mypyArgs": ["--follow-imports=normal", "--show-column-numbers"],
        "python.testing.pytestEnabled": true,
        "python.testing.pytestArgs": ["-o", "junit_family=xunit1"],
    }


Set up the documentation build as a task in ``.vscode/tasks.json``:

.. code-block:: json

    {
        "version": "2.0.0",
        "tasks": [
            {
                "label": "Build document",
                "type": "process",
                "options": {
                    "cwd": "${workspaceFolder}/docs"
                },
                "command": "${config:python.pythonPath}",
                "args": ["-m", "sphinx", "-b", "html", ".", "_build/html"],
                "group": {
                    "kind": "build",
                    "isDefault": true
                },
                "presentation": {
                    "echo": true,
                    "reveal": "silent",
                    "focus": false,
                    "panel": "dedicated",
                    "showReuseMessage": true,
                    "clear": true
                }
            }
        ]
    }

Set up the the debug shortcut in ``.vscode/launch.json``:

.. code-block:: json

    {
        "version": "0.2.0",
        "configurations": [
            {
                "name": "Run CharGer",
                "type": "python",
                "request": "launch",
                "program": "${workspaceFolder}/scripts/debug_example.py"
            }
        ]
    }

.. _Visual Studio Code: https://code.visualstudio.com/
