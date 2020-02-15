from typing import List

from loguru import logger
from typing_extensions import Final

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

        >>> config = CharGerConfig(...)
        >>> charger = CharGer()
        >>> charger.setup()
    """

    def __init__(self, config: CharGerConfig):
        self.config: Final[CharGerConfig] = config
        """Configuration as a :py:class:`~charger.config.CharGerConfig` object."""
        self.input_variants: List[Variant] = []
        """Parsed input variants."""

    def setup(self) -> None:
        """Setup all the intput data and annotation.

        Sequentially it calls:

            1. :py:meth:`_read_input_vcf`
        """
        self._validate_config()
        self._read_input_vcf()
        self._read_inheritance_gene_list()

    def _validate_config(self) -> None:
        """Validate the configuration.

        This method is automatically called when CharGer is created."""
        logger.info(f"Validate the given config")
        logger.debug(f"Given config: {self.config!r}")

    def _read_input_vcf(self) -> None:
        """Read input VCF.

        Input VCF is read from :py:attr:`.CharGerConfig.input`.
        """
        if self.config.input is None:
            raise ValueError(f"No input file is given in the config")

        logger.info(f"Read input VCF from {self.config.input}")
        self.input_variants = list(Variant.read_vcf(self.config.input, parse_csq=True))
        logger.info(
            f"Read total {len(self.input_variants)} variants from the input VCF"
        )

    def _read_inheritance_gene_list(self):
        """Read inheritance gene list"""
        if self.config.inheritance_gene_list is not None:
            self.config.inheritance_gene_list
