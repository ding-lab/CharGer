"""
CharGer (Characterization of Germline variants) is a software tool for interpreting and predicting clinical pathogenicity of germline variants.
"""  # noqa
import argparse
import sys

from loguru import logger

from .argtype import PathType
from .config import CharGerConfig

logger.disable("charger")

# Additional console help message text at the end
epilog = """
If you use CharGer, please cite our publication so we can continue to support CharGer development: 

    Scott, A.D. et al. (2018). Bioinformatics. https://doi.org/10.1093/bioinformatics/bty649.
"""  # noqa


def create_console_parser() -> argparse.ArgumentParser:
    """Create the CLI parser."""

    class ConsoleHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
    ):
        pass

    # Obtain the config defaults
    defaults = CharGerConfig()

    parser = argparse.ArgumentParser(
        description=__doc__, epilog=epilog, formatter_class=ConsoleHelpFormatter
    )
    parser.add_argument(
        "--input",
        metavar="VCF",
        type=PathType(exists=True),
        help="Path to the input VCF to be annotated",
    )
    parser.add_argument(
        "--output", metavar="TSV", type=PathType(), help="Path to CharGer output",
    )

    parser.add_argument(
        "--disease-specific",
        action="store_true",
        help="Enable disease specific inheritance-gene-list detection",
    )
    parser.add_argument(
        "--inheritance-gene-list",
        type=PathType(exists=True),
        metavar="TSV",
        help=(
            "Path to inheritance gene tab separated table. "
            "Table columns must be: gene, disease, mode_of_inheritance"
        ),
    )
    parser.add_argument(
        "--pathogenic-variant",
        metavar="VCF",
        type=PathType(exists=True),
        help="Path to the known pathogenic variants",
    )
    parser.add_argument(
        "--hotspot3d-cluster",
        metavar="TSV",
        type=PathType(exists=True),
        help="Path to HotSpot3D clusters result",
    )
    parser.add_argument(
        "--override-variant-info",
        action="store_true",
        help="Override the variant info using ClinVar description",
    )
    parser.add_argument(
        "--include-vcf-details",
        action="store_true",
        help="Include the VCF details in the output",
    )
    parser.add_argument(
        "--PP2-gene-list",
        type=PathType(exists=True),
        metavar="TXT",
        help="Path to PP2 gene list (list of gene symbols)",
    )

    anno_src_grp = parser.add_argument_group("annotation sources")
    anno_src_grp.add_argument(
        "--use-clinvar", action="store_true", help="Use ClinVar",
    )
    anno_src_grp.add_argument(
        "--clinvar-src",
        metavar="TSV",
        type=PathType(exists=True),
        help="Path to the ClinVar database to look up input variants",
    )

    threshold_grp = parser.add_argument_group("thresholds")
    threshold_grp.add_argument(
        "--rare-threshold",
        type=float,
        metavar="FREQ",
        default=defaults.rare_threshold,
        help="Maximal allele frequency to be a rare variant",
    )
    threshold_grp.add_argument(
        "--common-threshold",
        type=float,
        metavar="FREQ",
        default=defaults.common_threshold,
        help="Minimal allele frequency to be a common variant",
    )

    cg_cls_threshold_grp = parser.add_argument_group(
        "CharGer variant classiciation thresholds"
    )
    cg_cls_threshold_grp.add_argument(
        "--min-pathogenic-score",
        type=int,
        metavar="INT",
        default=defaults.min_pathogenic_score,
        help="Minimal total score for a variant to be pathogenic",
    )
    cg_cls_threshold_grp.add_argument(
        "--min-likely-pathogenic-score",
        type=int,
        metavar="INT",
        default=defaults.min_likely_pathogenic_score,
        help="Minimal total score for a variant to be likely pathogenic",
    )
    cg_cls_threshold_grp.add_argument(
        "--max-likely-benign-score",
        type=int,
        metavar="INT",
        default=defaults.max_likely_benign_score,
        help="Maximal total score for a variant to be likely benign",
    )
    cg_cls_threshold_grp.add_argument(
        "--max-benign-score",
        type=int,
        metavar="INT",
        default=defaults.max_benign_score,
        help="Maximal total score for a variant to be benign",
    )
    return parser


def parse_console() -> CharGerConfig:
    parser = create_console_parser()
    config = parser.parse_args(namespace=CharGerConfig())
    return config


def run() -> None:
    # Set up stderr format
    logger.remove()
    logger.add(
        sys.stderr,
        format=(
            # "<green>{time:YYYY-MM-DD HH:mm:ss}</green> "
            "<level>{level: <8}</level> "
            "| <level>{message}</level>"
        ),
    )
    # By default all the logging messages are disabled
    logger.enable("charger")
    config = parse_console()
    logger.info(f"Current config: {config!r}")
