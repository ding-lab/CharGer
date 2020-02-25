## Installation
Make sure the conda is available and set up [bioconda]. Create a new conda environment (for example, `charger_py37`):

    conda create -n charger_py37 python=3.7 pip poetry cyvcf2 pysam
    conda activate charger_py37
    git clone https://github.com/ding-lab/LightCharGer
    poetry install      # install charger and all its dependencies
    pre-commit install  # enforce style check at every commit

For python3.8, one needs to build the cyvcf2 manually from source:

    conda create -n charger_py38 python=3.8 poetry pip numpy
    conda activate charger_py38
    git clone https://github.com/ding-lab/LightCharGer
    # Install cyvcf2 from source, can skip once bioconda builds python3.8 cyvcf
    # conda create -n cyvcf python=3.8 cython numpy wheel
    # conda activate cyvcf
    # git clone --recursive https://github.com/brentp/cyvcf2
    # cd cyvcf2/htslib
    # autoheader
    # autoconf
    # ./configure --enable-libcurl
    # make
    # cd ..
    # python setup.py bdist_wheel
    # Additionally, run python setup.py test to make sure it passes all the tests
    # a .whl cyvcf2 package will be available under dist/; can be re-used
    pip install dist/cyvcf2-0.11.6-cp38-cp38-macosx_10_9_x86_64.whl
    poetry install
    pre-commit install

[bioconda]: https://bioconda.github.io/

## Usage

    charger -h

Visit CharGer's documentations for its detailed usage.


## Development

### Pass tests before commit
Style checks are enforced to pass before any git commit (by running `pre-commit install` once).
To run the style checks at any time, run the following command:

    pre-commit run -a

Otherwise, run the style checks manually by:

    isort               # Sort the import order
    black src tests     # Format the code
    flake8              # Check coding style

Functional tests and type checks are not run automatically for every commit,
but it's recommended to run and pass all the tests and type checks by:

    pytest -v
    mypy --pretty src tests

The repo should always pass all the tests described above.

### Build documentation
CharGer's documentation is powered by [sphinx] under `docs`. Build a new version by:

    cd docs
    make html

And the documentation will be available under `docs/_build/html`.


### Developing using Visual Studio Code
Here are some additional setup that utilize [Visual Studio Code]'s IDE:

- Run style checks ([black], [mypy] and [flake8]) at every file save
- The default build task will call [sphinx] to build the documentation
- A debug shortcut to go into CharGer internals using a startup script `scripts/debug_example.py`

Add the following workspace settings `.vscode/settings.json`:

```json
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
```

Set up the documentation build as a task in `.vscode/tasks.json`:

```json
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
```

Set up the the debug shortcut in `.vscode/launch.json`:

```json
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
```

[black]: https://github.com/psf/black
[mypy]: http://www.mypy-lang.org/
[flake8]: https://flake8.pycqa.org/
[sphinx]: https://www.sphinx-doc.org/
[Visual Studio Code]: https://code.visualstudio.com/
