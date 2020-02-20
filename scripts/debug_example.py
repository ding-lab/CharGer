from pathlib import Path

from charger.classifier import CharGer
from charger.config import CharGerConfig
from charger.console import setup_logger

REPO_ROOT = Path(__file__).parent.parent

setup_logger()
config = CharGerConfig(
    input=REPO_ROOT.joinpath("tests/examples/grch37_vep85_5_variants.vcf"),
    pathogenic_variant=REPO_ROOT.joinpath(
        "tests/examples/annotations/grch37_pathogenic_variants.vcf.gz"
    ),
    inheritance_gene_table=REPO_ROOT.joinpath(
        "tests/examples/annotations/inheritance_gene_table.tsv.gz"
    ),
    PP2_gene_list=REPO_ROOT.joinpath("tests/examples/annotations/pp2_gene_list.txt.gz"),
)
charger = CharGer(config)
charger.setup()
