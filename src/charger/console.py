import argparse
from pathlib import Path


class CharGerConfig(argparse.Namespace):
    """CharGer configuration."""

    # Define all the config options and their types

    @classmethod
    def create_console_parser(cls) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        return parser


def run() -> None:
    pass
