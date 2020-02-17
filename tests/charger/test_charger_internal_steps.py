import pytest

from charger.charger import CharGer
from charger.config import CharGerConfig
from charger.variant import VariantInheritanceMode


@pytest.fixture
def grch38(test_root):
    config = CharGerConfig(
        input=test_root.joinpath("examples/grch38_vep95_50_variants.vcf.gz"),
        inheritance_gene_list=test_root.joinpath(
            "examples/annotations/inheritance_gene_list.tsv.gz"
        ),
    )
    return CharGer(config)


def test_read_input(grch38):
    grch38._read_input_vcf()
    assert len(grch38.input_variants) == 50


def test_read_inheritance_gene_list(grch38):
    grch38._read_inheritance_gene_list()
    assert len(grch38.inheritance_genes) == 152
    assert grch38.inheritance_genes["BRCA1"] == VariantInheritanceMode.AUTO_DOMINANT
    assert (
        grch38.inheritance_genes["TERT"]
        == VariantInheritanceMode.AUTO_DOMINANT | VariantInheritanceMode.AUTO_RECESSIVE
    )


def test_read_inheritance_gene_list_disease_specific(test_root):
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
