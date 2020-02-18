from enum import Enum, auto
from typing import Any, Callable, Dict, List, Optional, Set

import attr
from loguru import logger
from typing_extensions import Final

from .config import ACMG_MODULES, CHARGER_MODULES, CharGerConfig
from .io import read_lines, read_tsv
from .variant import Variant, VariantInheritanceMode

logger.disable("charger")  # Disable emit logs by default


class CharGer:
    """
    Variant classifier.

    Args:
        config: CharGer's configurations.
            See :class:`~charger.config.CharGerConfig` for details to set it up.

    Examples:

        >>> config = CharGerConfig(...)
        >>> charger = CharGer()
        >>> charger.setup()
    """

    def __init__(self, config: CharGerConfig):
        self.config: Final[CharGerConfig] = config
        """Configuration as a :class:`~charger.config.CharGerConfig` object."""

        self.input_variants: List[Variant] = []
        """Parsed input variants."""

        self.pathogenic_variants: List[Variant] = []
        """Known pathogenic variants."""

        self.inheritance_genes: Dict[str, Optional[VariantInheritanceMode]] = {}
        """Variant inheritance dominance mode of the genes for PVS1 module."""

        self.pp2_genes: Set[str] = set()
        """Genes marked for PP2 module."""

        self.bp1_genes: Set[str] = set()
        """Genes marked for BP1 module."""

        self.results: List[CharGerResult] = []
        """
        Classification results of the input variants.
        The length and order should always be the same as :attr:`input_variants`.
        """

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
        """The availability of all CharGer modules. By default all modules are active."""

    def setup(self) -> None:
        """Setup all the intput data and annotation.

        Sequentially it calls:

            1. :meth:`_validate_config`
            2. :meth:`_read_input_variants`
            3. :meth:`_read_pathogenic_variants`
            4. :meth:`_read_inheritance_gene_list`
            5. :meth:`_read_pp2_gene_list`
            6. :meth:`_read_bp1_gene_list`
        """
        self._validate_config()
        self._read_input_variants()
        self._read_pathogenic_variants()
        self._read_inheritance_gene_list()
        self._read_pp2_gene_list()
        self._read_bp1_gene_list()

    def _validate_config(self) -> None:
        """Validate the configuration."""
        logger.info(f"Validate the given config")
        logger.debug(f"Given config: {self.config!r}")

    def _read_input_variants(self) -> None:
        """Read input VCF and set up the result template

        Load :attr:`input_variants` from :attr:`self.config.input <.CharGerConfig.input>`.
        Also populate :attr:`results` matching the input variant.
        """
        if self.config.input is None:
            raise ValueError(f"No input file is given in the config")

        logger.info(f"Read input VCF from {self.config.input}")
        # TODO: Skip variants with filter, or with high allele frequency
        # num_skipped_variants: Dict[str, int] = {"has_filter": 0}
        for variant in Variant.read_vcf(self.config.input, parse_csq=True):
            # # Skip the variant with filter (not PASS)
            # if variant.filter:
            #     logger.warning(
            #         f"{variant} has filter {','.join(variant.filter)}. Skipped"
            #     )
            #     num_skipped_variants["has_filter"] += 1
            #     continue
            self.input_variants.append(variant)

            # We also create the result template
            self.results.append(CharGerResult(variant))

        logger.info(
            f"Read total {len(self.input_variants)} variants from the input VCF"
        )

    def _read_pathogenic_variants(self) -> None:
        """Read known pathogenic variants.

        Load :attr:`pathogenic_variants`
        from :attr:`self.config.pathogenic_variant <.CharGerConfig.pathogenic_variant>`.
        """
        if self.config.pathogenic_variant is None:
            return
        logger.info(f"Read pathogenic VCF from {self.config.pathogenic_variant}")
        self.pathogenic_variants = list(
            Variant.read_vcf(self.config.pathogenic_variant, parse_csq=True)
        )
        logger.info(
            f"Read total {len(self.pathogenic_variants)} pathogenic variants from the VCF"
        )

    def _read_inheritance_gene_list(self) -> None:
        """Read inheritance gene list for PVS1 module.

        Load :attr:`inheritance_genes`
        from :attr:`self.config.inheritance_gene_list <.CharGerConfig.inheritance_gene_list>`.
        Skip PVS1 module if it's not provided.
        """
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

    def _read_pp2_gene_list(self) -> None:
        """Read gene list for PP2 module.

        Load :attr:`pp2_genes`
        from :attr:`self.config.PP2_gene_list <.CharGerConfig.PP2_gene_list>`.
        Skip PP2 module if not provided.
        """
        gene_list_pth = self.config.PP2_gene_list
        # Disable PP2 module if no list is provided
        if gene_list_pth is None:
            logger.warning(
                "CharGer cannot make PP2 calls without the given gene list. "
                "Disable PP2 module"
            )
            self._acmg_module_availability["PP2"] = ModuleAvailability.INVALID_SETUP
            return

        logger.info(f"Read PP2 gene list from {gene_list_pth}")
        self.pp2_genes = set(l.strip() for l in read_lines(gene_list_pth))
        logger.info(f"Marked {len(self.pp2_genes)} genes for PP2")

    def _read_bp1_gene_list(self) -> None:
        """Read gene list for BP1 module.

        Load :attr:`bp1_genes`
        from :attr:`self.config.BP1_gene_list <.CharGerConfig.BP1_gene_list>`.
        Skip BP1 module if not provided.
        """
        gene_list_pth = self.config.BP1_gene_list
        # Disable BP1 module if no list is provided
        if gene_list_pth is None:
            logger.warning(
                "CharGer cannot make BP1 calls without the given gene list. "
                "Disable BP1 module"
            )
            self._acmg_module_availability["BP1"] = ModuleAvailability.INVALID_SETUP
            return

        logger.info(f"Read BP1 gene list from {gene_list_pth}")
        self.bp1_genes = set(l.strip() for l in read_lines(gene_list_pth))
        logger.info(f"Marked {len(self.bp1_genes)} genes for BP1")


class ModuleAvailability(Enum):
    """Availability of a  ACMG or CharGer modules.

    Used by :attr:`CharGer._acmg_module_availability` and :attr:`CharGer._charger_module_availability`.

    Examples:

        Skip PVS1 module by:

        >>> charger = CharGer(CharGerConfig())
        >>> charger._acmg_module_availability = ModuleAvailability.USER_DISABLED
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


class ModuleDecision(Enum):
    """The decision of a ACMG or CharGer module on one variant.

    Used by :attr:`CharGerResult.acmg_decisions` and :attr:`CharGerResult.charger_decisions`
    """

    PASSED = auto()
    """The variant matched the criteria of the module."""
    FAILED = auto()
    """The variant didn't match the criteria of the module."""
    SKIPPED = auto()
    """The module was skipped."""

    @classmethod
    def _gen_decision_template(
        cls, available_modules: Dict[str, List[str]]
    ) -> Callable[[], Dict[str, Optional["ModuleDecision"]]]:
        """Generate the decision template for all the available modules."""

        def gen_template():
            decisions = {}
            for module_type, modules in available_modules.items():
                for m in modules:
                    decisions[m] = None
            return decisions

        return gen_template


@attr.s(auto_attribs=True, eq=False, order=False, slots=True)
class CharGerResult:
    """Result of the CharGer classification."""

    variant: Variant
    """The input variant."""

    clinvar: Dict[str, Any] = attr.Factory(dict)
    """Clinvar annotation."""

    acmg_decisions: Dict[str, Optional[ModuleDecision]] = attr.Factory(
        ModuleDecision._gen_decision_template(ACMG_MODULES)
    )
    """
    The decision of each ACMG module of the variant.

    `None` if the module is not run. See :class:`ModuleDecision` for the possible decisions.
    """

    charger_decisions: Dict[str, Optional[ModuleDecision]] = attr.Factory(
        ModuleDecision._gen_decision_template(CHARGER_MODULES)
    )
    """
    The decision of each CharGer module of the variant.

    Same usage as :attr:`acmg_decisions`.
    """
