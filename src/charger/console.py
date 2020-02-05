"""
CharGer (Characterization of Germline variants) is a software tool for interpreting and predicting clinical pathogenicity of germline variants.
"""  # noqa
import argparse
from pathlib import Path
import sys
from loguru import logger
from .argtype import PathType

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
        default=0.0005,
        help="Maximal allele frequency to be a rare variant",
    )
    threshold_grp.add_argument(
        "--common-threshold",
        type=float,
        metavar="FREQ",
        default=0.005,
        help="Minimal allele frequency to be a common variant",
    )
    return parser


class CharGerConfig(argparse.Namespace):
    """CharGer configuration."""

    # Define all the config options and their types
    input: Path
    output: Path
    disease_specific: bool
    inheritance_gene_list: Path
    pathogenic_variant: Path
    override_variant_info: bool
    include_vcf_details: bool
    PP2_gene_list: Path
    # annotation sources:
    use_clinvar: bool
    clinvar_src: Path
    rare_threshold: float
    common_threshold: float

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _get_kwargs(self):
        # Don't sort the arguments alphabetically
        return self.__dict__.items()

    def __repr__(self):
        # Modified from argparse._AttributeHolder.__repr__
        type_name = type(self).__name__
        arg_strings = []
        star_args = {}
        for arg in self._get_args():
            arg_strings.append(f"    {repr(arg)}")
        for name, value in self._get_kwargs():
            if name.isidentifier():
                arg_strings.append(f"    {name!s}={value!r}")
            else:
                star_args[name] = value
        if star_args:
            arg_strings.append(f"  **{repr(star_args)}")
        return "%s(\n%s\n)" % (type_name, ",\n".join(arg_strings))


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
    parser = create_console_parser()
    config = parser.parse_args(namespace=CharGerConfig())
    logger.info(f"Charger config:\n{config!r}")
