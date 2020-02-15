import logging
from pathlib import Path

import pytest
from _pytest.logging import caplog as _caplog  # noqa: F401
from loguru import logger

from charger.config import CharGerConfig


@pytest.fixture(scope="session")
def test_root():
    return Path(__file__).parent


@pytest.fixture
def caplog(_caplog):  # noqa: F811
    """A fixture to capture loguru logging messages.

    Copied from https://loguru.readthedocs.io/en/stable/resources/migration.html
    """

    class PropagateHandler(logging.Handler):
        def emit(self, record):
            logging.getLogger(record.name).handle(record)

    logger.enable("charger")
    handler_id = logger.add(PropagateHandler(), format="{message}")
    yield _caplog
    logger.remove(handler_id)
    logger.disable("charger")


@pytest.fixture(scope="session")
def default_config() -> CharGerConfig:
    """Default CharGer config."""
    return CharGerConfig()


@pytest.fixture(scope="session")
def example_input_vcf(test_root) -> Path:
    return test_root.joinpath("examples/grch38_vep95_50_variants.vcf.gz")
