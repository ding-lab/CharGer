## Installation

    conda create -n charger_py37 python=3.7 pip poetry pre-commit cyvcf2
    conda activate charger_py37
    poetry install
    pre-commit install  # Enable style check


For python3.8, one needs to build the cyvcf2 manually from source:

    conda create -n charger_py38 python=3.8 poetry pip cython
    conda activate charger_py38
    # Install cyvcf2 from source
    poetry install


## Usage

    charger -h

## Development
Run all the tests and type checks by:

    pytest -v
    mypy --pretty src

Style checks are enforced before any git commit using pre-commit. Run the style checks at any time by:

    pre-commit -a

Otherwise, run the style checks manually by:

    isort               # Sort the import order
    black src tests     # Format the code
    flake8              # Check coding style

The repo should always pass all the tests described above.

### Build documentation
CharGer's documentation is powered by sphinx under `docs`. Build a new version by:

    cd docs
    make html

And the documentation will be available under `docs/_build/html`.


### Developing using VS Code
The following VSCode workspace settings will run black, mypy and flake8 at every file save:

```json
{
    "python.formatting.provider": "black",
    "editor.formatOnSave": true,
    "editor.wordWrapColumn": 120,
    "python.linting.enabled": true,
    "python.linting.flake8Enabled": true,
    "python.linting.mypyEnabled": true,
    "python.testing.pytestEnabled": true,
    "python.testing.pytestArgs": ["-o", "junit_family=xunit1"],
}
```


