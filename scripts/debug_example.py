from pathlib import Path

from charger.charger import CharGer
from charger.config import CharGerConfig
from charger.console import setup_logger

REPO_ROOT = Path(__file__).parent.parent

setup_logger()
config = CharGerConfig(
    input=REPO_ROOT.joinpath("tests/examples/grch38_vep95_50_variants.vcf.gz"),
    inheritance_gene_list=REPO_ROOT.joinpath(
        "tests/examples/annotations/inheritance_gene_list.tsv.gz"
    ),
    PP2_gene_list=REPO_ROOT.joinpath("tests/examples/annotations/pp2_gene_list.txt.gz"),
)
charger = CharGer(config)
charger.setup()
