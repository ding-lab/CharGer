from enum import Enum, auto
from typing import Dict, List, Optional

from loguru import logger
from typing_extensions import Final

from .config import ACMG_MODULES, CHARGER_MODULES, CharGerConfig
from .io import read_tsv
from .variant import Variant, VariantInheritanceMode

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

        self.inheritance_genes: Dict[str, Optional[VariantInheritanceMode]] = {}
        """Given variant inheritance dominance mode of the genes."""

        self._acmg_module_availability: Dict[str, ModuleAvailability] = {
            m: ModuleAvailability.ACTIVE
            for module_type, modules in ACMG_MODULES.items()
            for m in modules
        }
        """The availability of all ACMG modules. By default all modules are active."""

        self._charger_module_availability: Dict[str, ModuleAvailability] = {
            m: ModuleAvailability.ACTIVE
            for module_type, modules in CHARGER_MODULES.items()
            for m in modules
        }
        """The availability of all CharGer modules. By default all modules are active"""

    def setup(self) -> None:
        """Setup all the intput data and annotation.

        Sequentially it calls:

            1. :py:meth:`_read_input_vcf`
            2. :py:meth`_read_inheritance_gene_list`
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

    def _read_inheritance_gene_list(self) -> None:
        """Read inheritance gene list"""
        tsv_pth = self.config.inheritance_gene_list
        # Disable PVS1 module if no list is provided
        if tsv_pth is None:
            logger.warning(
                "CharGer cannot make PVS1 calls without inheritance gene list. "
                "Disable PVS1 module"
            )
            self._acmg_module_availability["PVS1"] = ModuleAvailability.INVALID_SETUP
            return
        if self.config.disease_specific:
            raise NotImplementedError(
                "Cannot read disease specific inheritance gene list"
            )

        logger.info(f"Read inheritance gene list from {tsv_pth}")
        reader = read_tsv(tsv_pth, as_dict=False)
        header = next(reader)
        if len(header) < 3:
            logger.error(
                f"Expect inheritance gene list to have at least three columns; "
                f"only got {', '.join(header)}"
            )
            raise ValueError(f"Invalid table format in {tsv_pth}")
        for (gene, diseases, modes_of_inheritance, *_) in reader:
            modes = VariantInheritanceMode.parse(modes_of_inheritance)
            self.inheritance_genes[gene] = modes
        logger.info(f"Loaded inheritance mode of {len(self.inheritance_genes)} genes")


class ModuleAvailability(Enum):
    """Availability of the ACMG and CharGer modules.

    Used by :py:attr:`CharGer._acmg_module_availability` and :py:attr:`CharGer._charger_module_availability`
    """

    ACTIVE = auto()
    """The module is active."""
    USER_DISABLED = auto()
    """The module is disabled by user. The module will be skipped."""
    INVALID_SETUP = auto()
    """
    The module has invalid setup, such as no additional annotation provided.
    The module will be skipped.
    """
