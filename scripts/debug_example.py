from pathlib import Path

from charger.classifier import CharGer
from charger.config import CharGerConfig
from charger.console import setup_logger

REPO_ROOT = Path(__file__).parent.parent
path_in_repo = REPO_ROOT.joinpath

# Enable logging to stderr
setup_logger()

# Configure CharGer
config = CharGerConfig(
    input=path_in_repo(
        "tests/examples/10.1056_NEJMoa1508054_S4_AD_vep85.sorted.vcf.gz"
    ),
    pathogenic_variant=path_in_repo(
        "tests/examples/annotations/grch37_pathogenic_variants.vcf.gz"
    ),
    inheritance_gene_table=path_in_repo(
        "tests/examples/annotations/inheritance_gene_table.tsv.gz"
    ),
    PP2_gene_list=path_in_repo("tests/examples/annotations/pp2_gene_list.txt.gz"),
    clinvar_table=path_in_repo(
        "tests/examples/annotations/clinvar_chrom_22_only.b37.tsv.gz"
    ),
)

# Initiate and run the CharGer classifier
charger = CharGer(config)
charger.setup()
charger.match_clinvar()
charger.run_acmg_modules()
charger.run_charger_modules()
