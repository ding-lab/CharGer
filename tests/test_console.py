from charger.console import CharGerConfig


def test_direct_config_creation():
    config = CharGerConfig(output=None)
    assert config.output is None
