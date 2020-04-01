from enum import Enum, auto
from typing import Any, Dict, List, Optional, Set

from loguru import logger
from pysam import TabixFile

from .acmg_modules import run_pm4, run_pvs1
from .config import ACMG_MODULES, CHARGER_MODULES, CharGerConfig
from .custom_modules import run_pmc1, run_ppc1, run_ppc2, run_psc1
from .io import read_lines, read_tsv
from .result import CharGerResult
from .variant import ClinicalSignificance, GeneInheritanceMode, Variant

try:
    from typing import Final
except ImportError:
    # Backport additional typings prior to python 3.8
    from typing_extensions import Final  # type: ignore


logger.disable("charger")  # Disable emit logs by default


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


InheritanceGenesType = Dict[str, Optional[GeneInheritanceMode]]
"""Type hint alias for :attr:`CharGer.inheritance_genes`."""


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
        >>> charger.match_clinvar()
        >>> charger.run_acmg_modules()
        >>> charger.run_charger_modules()
    """

    def __init__(self, config: CharGerConfig):
        self.config: Final[CharGerConfig] = config
        """Configuration as a :class:`~charger.config.CharGerConfig` object."""

        self.input_variants: List[Variant] = []
        """Parsed input variants."""

        self.pathogenic_variants: List[Variant] = []
        """Known pathogenic variants."""

        self.inheritance_genes: InheritanceGenesType = {}
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
            4. :meth:`_read_inheritance_gene_table`
            5. :meth:`_read_pp2_gene_list`
            6. :meth:`_read_bp1_gene_list`
        """
        self._validate_config()
        self._read_input_variants()
        self._read_pathogenic_variants()
        self._read_inheritance_gene_table()
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
        for variant in Variant.read_and_parse_vcf(self.config.input):
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

        logger.success(
            f"Read total {len(self.input_variants):,d} variants from the input VCF"
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
            Variant.read_and_parse_vcf(self.config.pathogenic_variant)
        )
        logger.info(
            f"Read total {len(self.pathogenic_variants):,d} pathogenic variants from the VCF"
        )

    def _read_inheritance_gene_table(self) -> None:
        """Read inheritance gene table for PVS1 module.

        Load :attr:`inheritance_genes`
        from :attr:`self.config.inheritance_gene_table <.CharGerConfig.inheritance_gene_table>`.
        Skip PVS1 module if it's not provided.
        """
        tsv_pth = self.config.inheritance_gene_table
        # Disable PVS1 module if no list is provided
        if tsv_pth is None:
            logger.warning(
                "Inheritance gene table is not provided, CharGer cannot make ACMG PVS1/PM4 calls "
                "or CharGer PSC1/PPC1 calls. Disable all these modules"
            )
            for m in ["PVS1", "PM4"]:
                self._acmg_module_availability[m] = ModuleAvailability.INVALID_SETUP
            for m in ["PSC1", "PPC1"]:
                self._charger_module_availability[m] = ModuleAvailability.INVALID_SETUP
            return
        else:
            logger.info(
                "Disable CharGer PMC1/PPC2 modules when inheritance gene table is provided"
            )
            for m in ["PMC1", "PPC2"]:
                self._charger_module_availability[m] = ModuleAvailability.INVALID_SETUP
        if self.config.disease_specific:
            raise NotImplementedError(
                "Cannot read disease specific inheritance gene table"
            )

        logger.info(f"Read inheritance gene table from {tsv_pth}")
        reader = read_tsv(tsv_pth, as_dict=False)
        header = next(reader)
        if len(header) < 3:
            logger.error(
                f"Expect inheritance gene table to have at least three columns; "
                f"only got {', '.join(header)}"
            )
            raise ValueError(f"Invalid table format in {tsv_pth}")
        for (gene, diseases, modes_of_inheritance, *_) in reader:
            modes = GeneInheritanceMode.parse(modes_of_inheritance)
            self.inheritance_genes[gene] = modes
        logger.info(
            f"Loaded inheritance mode of {len(self.inheritance_genes):,d} genes"
        )

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
        logger.info(f"Marked {len(self.pp2_genes):,d} genes for PP2")

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
        logger.info(f"Marked {len(self.bp1_genes):,d} genes for BP1")

    @staticmethod
    def _match_clinvar_one_variant(
        variant: Variant, tabix: TabixFile, cols: List[str]
    ) -> Optional[Dict[str, Any]]:
        """Match the variant to the given ClinVar tabix table.

        Args:
            variant: Variant to be matched
            tabix: Tabix indexed CliVar table
            cols: All ClinVar columns in the table

        Returns:
            None if no ClinVar match. When matched, returns a `dict` of the clinvar record,
            where the key ``final_clinical_significance`` stores the final clinical significance type
            in :class:`ClinicalSignificance`.
        """
        try:
            # TabixFile.fetch will raise ValueError if the given region is out of bound
            row_iter = tabix.fetch(
                region=f"{variant.chrom}:{variant.start_pos}-{variant.end_pos}"
            )
        except ValueError as e:
            # Do nothing if it's querying for a chromosome not in the ClinVar table
            if "could not create iterator for region" not in e.args[0]:
                logger.opt(exception=e).debug(f"Tabix fetch ClinVar failed: {e}")
            return None

        for row in row_iter:
            record = dict(zip(cols, row.split("\t")))
            if (
                int(record["start"]) == variant.start_pos
                and int(record["stop"]) == variant.end_pos
                and record["alt"] == variant.alt_allele
            ):
                if record["ref"] != variant.ref_allele:
                    logger.warning(
                        f"{variant!r} got a clinvar match but their reference alleles are different: "
                        f"{variant.ref_allele!r} != {record['ref']!r}"
                    )
                # Parse the clinical significance of the record
                record[
                    "final_clinical_significance"
                ] = ClinicalSignificance.parse_clinvar_record(record)
                return record
        return None

    def match_clinvar(self) -> None:
        """Match the input variant with the ClinVar table.

        Update :attr:`CharGerResult.clinvar` the variant matches a ClinVar record
        by calling :meth:`_match_clinvar_one_variant`.
        """
        if self.config.clinvar_table is None:
            logger.info("Skip matching ClinVar")
            return
        logger.info(
            f"Match input variants with ClinVar table at {self.config.clinvar_table}"
        )
        clinvar_match_num = 0
        with TabixFile(str(self.config.clinvar_table), encoding="utf8") as tabix:
            cols = tabix.header[0][len("#") :].split("\t")
            for result in self.results:
                record = self._match_clinvar_one_variant(result.variant, tabix, cols)
                if record is not None:
                    result.clinvar = record
                    clinvar_match_num += 1
        logger.success(
            f"Matched {clinvar_match_num:,d} out of {len(self.input_variants):,d} input variants to a ClinVar record"
        )

    @staticmethod
    def _run_or_skip_module(
        module_name: str, module_avail: "ModuleAvailability"
    ) -> bool:
        if module_avail is ModuleAvailability.ACTIVE:
            logger.info("Running {name} module", name=module_name)
            return True
        else:
            logger.info("Skipped {name} module", name=module_name)
            return False

    def run_acmg_modules(self) -> None:
        """Run all ACMG modules.

        See :mod:`~charger.acmg_modules` for all the currently implemented
        modules.
        """
        logger.info("Run all ACMG modules")

        def run_or_skip(module_name: str):
            return self._run_or_skip_module(
                module_name, self._acmg_module_availability[module_name]
            )

        # PVS1
        if run_or_skip("PVS1"):
            for result in self.results:
                run_pvs1(result, self.inheritance_genes)

        # PM4
        if run_or_skip("PM4"):
            for result in self.results:
                run_pm4(result, self.inheritance_genes)

    def run_charger_modules(self) -> None:
        """Run all CharGer customized modules.

        See :mod:`~charger.custom_modules` for all the currently implemented
        modules.
        """
        logger.info("Run all CharGer modules")

        def run_or_skip(module_name: str):
            return self._run_or_skip_module(
                module_name, self._charger_module_availability[module_name]
            )

        # PSC1
        if run_or_skip("PSC1"):
            for result in self.results:
                run_psc1(result, self.inheritance_genes)

        # PMC1
        if run_or_skip("PMC1"):
            for result in self.results:
                run_pmc1(result, self.inheritance_genes)

        # PPC1
        if run_or_skip("PPC1"):
            for result in self.results:
                run_ppc1(result, self.inheritance_genes)

        # PPC
        if run_or_skip("PPC2"):
            for result in self.results:
                run_ppc2(result, self.inheritance_genes)
