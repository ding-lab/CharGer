from typing import List

from loguru import logger

from .config import CharGerConfig
from .variant import Variant

logger.disable("charger")  # Disable emit logs by default


class CharGer:
    """
    Variant classifier.

    Args:
        config: CharGer's configurations.
            See :py:class:`~charger.config.CharGerConfig` for details to set it up.

    Example:

        >>> charger = CharGer(CharGerConfig())
        >>> charger.read_input_data()
    """

    def __init__(self, config: CharGerConfig):
        self.config: CharGerConfig = config
        self.user_variants: List[Variant] = []
        self.validate_config()

    def validate_config(self):
        """Validate the configuration."""
        if not self.config.input.exists():
            raise ValueError(f"Input file does not exist at {self.config.input}")

    def read_input_data(self):
        """Read the input data.

        Input VCF is read from :py:attr:`.CharGerConfig.input`.
        """
        logger.debug(f"Read input VCF from {self.config.input}")
        self.user_variants = list(Variant.read_vcf(self.config.input, parse_csq=True))
        logger.info(f"Read total {len(self.user_variants)} variants from the input VCF")
