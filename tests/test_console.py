from charger.config import CharGerConfig
from charger.console import create_console_parser


def test_empty_config_creation():
    config = CharGerConfig()
    assert config.input is None
    assert config.output is None


def test_default_config():
    parser = create_console_parser()
    config = parser.parse_args([], namespace=CharGerConfig())
    assert config == CharGerConfig()
