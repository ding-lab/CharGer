from pathlib import Path

from charger.classifier import CharGer
from charger.config import CharGerConfig
from charger.console import setup_logger

REPO_ROOT = Path(__file__).parent.parent

setup_logger()
config = CharGerConfig(
    input=REPO_ROOT.joinpath(
        "tests/examples/10.1056_NEJMoa1508054_S4_AD_vep85.sorted.vcf.gz"
    ),
    pathogenic_variant=REPO_ROOT.joinpath(
        "tests/examples/annotations/grch37_pathogenic_variants.vcf.gz"
    ),
    inheritance_gene_table=REPO_ROOT.joinpath(
        "tests/examples/annotations/inheritance_gene_table.tsv.gz"
    ),
    PP2_gene_list=REPO_ROOT.joinpath("tests/examples/annotations/pp2_gene_list.txt.gz"),
    clinvar_table=REPO_ROOT.joinpath(
        "tests/examples/annotations/clinvar_chrom_22_only.b37.tsv.gz"
    ),
)
charger = CharGer(config)
charger.setup()
charger.match_clinvar()
charger.run_acmg_modules()
charger.run_charger_modules()
