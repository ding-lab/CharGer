from typing import List

import attr
import pytest

from charger.config import CharGerConfig
from charger.console import create_console_parser, parse_console


@pytest.fixture
def required_args(test_root):
    example_input_vcf = test_root.joinpath("examples/grch38_vep95_50_variants.vcf.gz")
    return [f"--input={example_input_vcf}"]


def test_default_config(default_config: CharGerConfig, required_args: List[str]):
    # Make sure we can launch CharGer with default settings
    parser = create_console_parser()
    config = parser.parse_args(required_args, namespace=CharGerConfig())

    # Remove the required arguments
    default_config_d = attr.asdict(default_config)
    config_d = attr.asdict(config)
    required_args = ["input"]
    for arg in required_args:
        default_config_d.pop(arg)
        config_d.pop(arg)
    assert config_d == default_config_d


def test_console_score_override(required_args: List[str], caplog):
    config = parse_console(
        required_args
        + [
            "--override-acmg-score=PS1=9 PS2=3",
            "--override-charger-score=PSC1=5 PMC1=3",
        ]
    )
    assert config.acmg_module_scores["PS1"] == 9
    assert config.acmg_module_scores["PS2"] == 3
    assert config.charger_module_scores["PSC1"] == 5
    assert config.charger_module_scores["PMC1"] == 3
    assert "Console parameters:" in caplog.text
    assert "--override-acmg-score=PS1=9 PS2=3" in caplog.text


def test_console_score_override_invalid_module(required_args: List[str]):
    with pytest.raises(SystemExit) as excinfo:
        parse_console(required_args + ["--override-acmg-score=PPAP=999"])
        assert "Module does not exist: PPAP" in excinfo.value.message
