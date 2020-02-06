from charger.console import CharGerConfig, create_console_parser


def test_empty_config_creation():
    CharGerConfig()
    assert CharGerConfig.input is None
    assert CharGerConfig.output is None


def test_default_config():
    parser = create_console_parser()
    config = parser.parse_args([], namespace=CharGerConfig())
    assert config == CharGerConfig()
