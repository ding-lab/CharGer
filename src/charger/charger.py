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
        >>> charger.read_input_variants()
    """

    def __init__(self, config: CharGerConfig):
        self.config: CharGerConfig = config
        self.user_variants: List[Variant] = []

    def read_input_variants(self):
        """Read the VEP-annotated variants from the input VCF.

        Input VCF is read from :py:attr:`.CharGerConfig.input`.
        """
