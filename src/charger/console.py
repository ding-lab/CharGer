"""
CharGer (Characterization of Germline variants) is a software tool for interpreting and predicting clinical pathogenicity of germline variants.
"""  # noqa
import argparse
from pathlib import Path
import sys
from typing import Optional
import attr
from loguru import logger
from .argtype import PathType

logger.disable("charger")

# Additional console help message text at the end
epilog = """
If you use CharGer, please cite our publication so we can continue to support CharGer development: 

    Scott, A.D. et al. (2018). Bioinformatics. https://doi.org/10.1093/bioinformatics/bty649.
"""  # noqa


@attr.s(auto_attribs=True, kw_only=True, repr=False)
class CharGerConfig:
    """CharGer configuration."""

    # Define all the config options and their types
    input: Optional[Path] = None
    output: Optional[Path] = None
    disease_specific: bool = False
    inheritance_gene_list: Optional[Path] = None
    hotspot3d_cluster: Optional[Path] = None
    pathogenic_variant: Optional[Path] = None
    override_variant_info: bool = False
    include_vcf_details: bool = False
    PP2_gene_list: Optional[Path] = None
    # annotation sources:
    use_clinvar: bool = False
    clinvar_src: Optional[Path] = None
    rare_threshold: float = 0.0005
    common_threshold: float = 0.005

    def __repr__(self):
        arg_strings = []
        for name, value in attr.asdict(self).items():
            arg_strings.append(f"    {name!s}={value!r},")
        return "CharGerConfig(\n{}\n)".format("\n".join(arg_strings))


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
