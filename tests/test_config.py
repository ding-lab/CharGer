from charger.config import CharGerConfig


def test_empty_config_creation() -> None:
    # Make sure we can always create an empty config
    config = CharGerConfig()
    assert config.input is None
    assert config.output is None
    assert config.acmg_module_scores["PVS1"] == 8
    assert config.charger_module_scores["BSC1"] == -6


def test_config_module_score_override(default_config: CharGerConfig) -> None:
    c = CharGerConfig()
    c.acmg_module_scores["PS1"] = 9
    c.charger_module_scores["BSC1"] = -9
    # Make sure each config has its own score dict
    assert c.acmg_module_scores is not default_config.acmg_module_scores
    assert c.charger_module_scores is not default_config.charger_module_scores
    assert c.acmg_module_scores["PS1"] != default_config.acmg_module_scores["PS1"]
    assert (
        c.charger_module_scores["BSC1"] != default_config.charger_module_scores["BSC1"]
    )
