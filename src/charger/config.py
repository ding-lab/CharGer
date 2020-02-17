from pathlib import Path
from typing import Dict, List, Optional, Tuple

import attr
from loguru import logger

logger.disable("charger")  # Disable emit logs by default


ACMG_MODULES: Dict[str, List[str]] = {
    # ACMG pathogenic modules
    "p_very_strong": ["PVS1"],
    "p_strong": ["PS1", "PS2", "PS3", "PS4"],
    "p_moderate": ["PM1", "PM2", "PM3", "PM4", "PM5", "PM6"],
    "p_support": ["PP1", "PP2", "PP3", "PP4", "PP5"],
    # ACMG benign modules
    "b_support": ["BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"],
    "b_moderate": [],
    "b_strong": ["BS1", "BS2", "BS3", "BS4"],
    "b_standalone": ["BA1"],
}

CHARGER_MODULES: Dict[str, List[str]] = {
    # CharGer pathogenic modules
    "p_strong": ["PSC1"],
    "p_moderate": ["PMC1"],
    "p_support": ["PPC1", "PPC2"],
    # CharGer benign modules
    "b_support": [],
    "b_moderate": ["BMC1"],
    "b_strong": ["BSC1"],
}


def _set_default_module_scores() -> Tuple[Dict[str, int], Dict[str, int]]:
    """Define the default scores for all ACMG and CharGer modules."""
    # Set up default module scores based on their type
    default_module_type_scores = {
        "p_very_strong": 8,
        "p_strong": 4,
        "p_moderate": 2,
        "p_support": 1,
        "b_support": -1,
        "b_moderate": -2,
        "b_strong": -4,
        "b_standalone": -8,
    }
    acmg_module_scores: Dict[str, int] = {}
    for module_type, modules in ACMG_MODULES.items():
        for m in modules:
            acmg_module_scores[m] = default_module_type_scores[module_type]

    charger_module_scores: Dict[str, int] = {}
    for module_type, modules in CHARGER_MODULES.items():
        for m in modules:
            charger_module_scores[m] = default_module_type_scores[module_type]
    # Override the score for special cases
    acmg_module_scores["PS1"] = 7  # ACMG p_very_strong
    charger_module_scores["BSC1"] = -6  # CHARGER b_strong
    return acmg_module_scores, charger_module_scores


# Create the default scores for char classification and each module
_default_acmg_scores, _default_charger_scores = _set_default_module_scores()


@attr.s(auto_attribs=True, kw_only=True, repr=True)
class CharGerConfig:
    """CharGer configuration.

    CharGer can be configured programatically using this object.

    Examples:

        >>> config = CharGerConfig()
    """

    # Define all the config options and their types
    input: Optional[Path] = None
    """Path to the input VCF."""

    output: Optional[Path] = None
    disease_specific: bool = False
    """Whether to enable disease specific detection for :attr:`inheritance_gene_list`."""

    inheritance_gene_list: Optional[Path] = None
    """Path to a TSV table of inheritance genes.

    The table columns must have these three columns in order: ``gene``, ``diseases``, ``modes_of_inheritance``.
    The first row is header and will be ignored.

        - gene: Gene symbol.
        - diseases: Name of the diseases seperated by space.
        - modes_of_inheritance: Name of the inheritance modes seperated by space. Possible modes:
          ``autosomal recessive``, ``autosomal dominant``,
          ``Y-linked``, ``X-linked recessive``, ``X-link dominant``,
          ``unknown``.

    .. csv-table:: Example table of ``inheritance_gene_list``
        :header-rows: 1
        :delim: |

        gene  | diseases     | modes_of_inheritance
        BRCA2 | BRCA, OV     | autosomal recessive, autosomal dominant
        NF1   | Neurofibroma | autosomal dominant
    """

    hotspot3d_cluster: Optional[Path] = None
    pathogenic_variant: Optional[Path] = None
    override_variant_info: bool = False
    include_vcf_details: bool = False

    # Annotation sources:
    use_clinvar: bool = False
    clinvar_src: Optional[Path] = None

    # Thresholds
    rare_threshold: float = 0.0005
    common_threshold: float = 0.005

    acmg_module_scores: Dict[str, int] = attr.Factory(_default_acmg_scores.copy)
    charger_module_scores: Dict[str, int] = attr.Factory(_default_charger_scores.copy)

    # CharGer variant classification thresholds
    min_pathogenic_score: int = 9
    min_likely_pathogenic_score: int = 5
    max_likely_benign_score: int = -4
    max_benign_score: int = -8

    # ACMG classifcation modules
    PP2_gene_list: Optional[Path] = None

    # CharGer classification modules
