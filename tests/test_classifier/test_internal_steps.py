from pathlib import Path

import pytest

from charger.classifier import CharGer, ModuleAvailability
from charger.config import CharGerConfig
from charger.variant import VariantInheritanceMode


@pytest.fixture
def grch38(test_root) -> CharGer:
    config = CharGerConfig(
        input=test_root.joinpath("examples/grch38_vep95_50_variants.vcf.gz"),
        pathogenic_variant=test_root.joinpath(
            "examples/annotations/grch37_pathogenic_variants.vcf.gz"
        ),
        inheritance_gene_list=test_root.joinpath(
            "examples/annotations/inheritance_gene_list.tsv.gz"
        ),
        PP2_gene_list=test_root.joinpath("examples/annotations/pp2_gene_list.txt.gz"),
    )
    return CharGer(config)


def test_read_input_variants(grch38: CharGer):
    grch38._read_input_variants()
    assert len(grch38.input_variants) == 50
    assert len(grch38.input_variants) == len(grch38.results)
    for v, result in zip(grch38.input_variants, grch38.results):
        assert v is result.variant


def test_read_pathogenic_variants(grch38: CharGer):
    grch38._read_pathogenic_variants()
    assert len(grch38.pathogenic_variants) == 1819


def test_read_inheritance_gene_list(grch38: CharGer):
    grch38._read_inheritance_gene_list()
    assert len(grch38.inheritance_genes) == 152
    assert grch38.inheritance_genes["BRCA1"] == VariantInheritanceMode.AUTO_DOMINANT
    assert (
        grch38.inheritance_genes["TERT"]
        == VariantInheritanceMode.AUTO_DOMINANT | VariantInheritanceMode.AUTO_RECESSIVE
    )


def test_no_inheritance_gene_list(grch38: CharGer):
    grch38.config.inheritance_gene_list = None
    grch38._read_inheritance_gene_list()
    assert grch38._acmg_module_availability["PVS1"] is ModuleAvailability.INVALID_SETUP


def test_read_inheritance_gene_list_disease_specific(test_root: Path):
    charger = CharGer(
        CharGerConfig(
            disease_specific=True,
            inheritance_gene_list=test_root.joinpath(
                "examples/annotations/inheritance_gene_list.tsv.gz"
            ),
        )
    )
    with pytest.raises(NotImplementedError):
        charger._read_inheritance_gene_list()


def test_read_pp2_gene_list(grch38: CharGer):
    grch38._read_pp2_gene_list()
    assert len(grch38.pp2_genes) == 152
    assert "RET" in grch38.pp2_genes
    assert "BRCA1" in grch38.pp2_genes


def test_no_bp1_gene_list(grch38: CharGer, caplog):
    grch38._read_bp1_gene_list()
    assert grch38._acmg_module_availability["BP1"] is ModuleAvailability.INVALID_SETUP
    assert "Disable BP1 module" in caplog.text
