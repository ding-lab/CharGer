import pytest

from charger.console import parse_console


def test_console_score_override():
    config = parse_console(
        ["--override-acmg-score=PS1=9 PS2=3", "--override-charger-score=PSC1=5 PMC1=3"]
    )
    assert config.acmg_module_scores["PS1"] == 9
    assert config.acmg_module_scores["PS2"] == 3
    assert config.charger_module_scores["PSC1"] == 5
    assert config.charger_module_scores["PMC1"] == 3


def test_console_score_override_invalid_module():
    with pytest.raises(SystemExit) as excinfo:
        parse_console(["--override-acmg-score=PPAP=999"])
        assert "Module does not exist: PPAP" in excinfo.value.message
