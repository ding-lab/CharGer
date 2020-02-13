from pathlib import Path

import pytest

from charger.variant import Variant

HERE = Path(__file__).parent


@pytest.fixture()
def snp(request):
    return Variant("1", 42, 42, "A", "G")


@pytest.fixture
def insertion():
    return Variant("1", 42, 42, "A", "ATAT")


@pytest.fixture
def deletion():
    return Variant("1", 42, 44, "ATA", "A")


@pytest.fixture
def sv():
    return Variant("1", 2827693, 2827701, "CCCCTCGCA", "C", info={"SVTYPE": "DEL"})


def test_variant_is_snp(snp, insertion, deletion, sv):
    assert snp.is_snp()
    assert not insertion.is_snp()
    assert not deletion.is_snp()
    assert not sv.is_snp()


def test_variant_is_sv(snp, insertion, deletion, sv):
    assert not snp.is_sv()
    assert not insertion.is_sv()
    assert not deletion.is_sv()
    assert sv.is_sv()


def test_variant_is_indel(snp, insertion, deletion, sv):
    assert not snp.is_indel()
    assert insertion.is_indel()
    assert deletion.is_indel()
    assert not sv.is_indel()


def test_variant_is_deletion(snp, insertion, deletion, sv):
    assert not snp.is_deletion()
    assert not insertion.is_deletion()
    assert deletion.is_deletion()
    assert not sv.is_deletion()


@pytest.fixture(scope="module")
def grch37_vep85_annotated_variants():
    return list(
        Variant.read_vcf(
            HERE / "examples" / "grch37_vep85_5_variants.vcf", parse_csq=True
        )
    )


@pytest.fixture(scope="module")
def grch38_vep95_annotated_variants():
    return list(
        Variant.read_vcf(
            HERE.joinpath("examples/grch38_vep95_50_variants.vcf.gz"), parse_csq=True,
        )
    )


@pytest.fixture(scope="module")
def grch38_vep95_annotated_variants_info_fixed():
    return list(
        Variant.read_vcf(
            HERE.joinpath("examples/grch38_vep95_50_variants.info_fixed.vcf.gz"),
            parse_csq=True,
        )
    )


@pytest.mark.parametrize(
    "vcf_pth,vep_ver",
    [
        ("examples/grch38_vep95_50_variants.bcf", 95),
        ("examples/grch38_vep95_50_variants.info_fixed.vcf.gz", 95),
        ("examples/grch38_vep95_50_variants.vcf.gz", 95),
        ("examples/grch37_vep85_5_variants.vcf", 85),
    ],
)
def test_read_vcf_detect_vep_version(caplog, vcf_pth, vep_ver):
    reader = Variant.read_vcf(HERE.joinpath(vcf_pth), parse_csq=True)
    next(reader)
    assert f"VEP version {vep_ver}" in caplog.text


def test_read_vcf_grch37(grch37_vep85_annotated_variants):
    variants = grch37_vep85_annotated_variants
    assert len(variants) == 5
    assert variants[0].chrom == "19"
    assert variants[0].start_pos == 45855804
    assert variants[0].end_pos == 45855805
    assert variants[0].ref_allele == "CT"
    assert variants[0].alt_allele == "C"
    assert variants[0].id is None
    assert variants[0].filter is None
    assert variants[0].is_indel()
    assert variants[0].is_deletion()
    assert not variants[0].is_snp()
    assert not variants[0].is_sv()


def test_read_vcf_grch38(grch38_vep95_annotated_variants):
    assert len(grch38_vep95_annotated_variants) == 50
    v = grch38_vep95_annotated_variants[42]
    assert v.info["AC"] == 2
    assert v.info["AF"] == 0.5
    assert v.info["AN"] == 4
    assert v.info["DP"] == 350
    assert v.info["set"] == "snp"
    assert v.info["CSQ"] is v.parsed_csq
    assert v.is_snp()
    assert not v.is_indel()
    assert not v.is_deletion()
    assert not v.is_sv()


def test_variant_parse_csq_grch38(grch38_vep95_annotated_variants):
    v = grch38_vep95_annotated_variants[42]
    assert len(v.parsed_csq) == 2
    csq0, csq1 = v.parsed_csq
    assert csq0["HGVSc"] == "ENST00000378785.6:c.414C>T"
    assert csq0["Consequence"] == "synonymous_variant"
    assert csq0["AF"] == "0.4251"
    assert csq0["ExAC_AF"] == ""
    assert csq0["gnomAD_AF"] == "0.3084"

    assert csq1["HGVSc"] == "ENST00000475091.2:c.258C>T"
    assert csq1["Consequence"] == "synonymous_variant"
    assert csq1["AF"] == "0.4251"
    assert csq1["ExAC_AF"] == ""
    assert csq1["gnomAD_AF"] == "0.3084"


def test_variant_parse_csq_grch37(grch37_vep85_annotated_variants):
    parsed_csq = grch37_vep85_annotated_variants[0].parsed_csq
    assert len(parsed_csq) == 11
    assert parsed_csq[10]["SYMBOL"] == "KLC3"
    assert parsed_csq[10]["ENSP"] == "ENSP00000466974"

    v = grch37_vep85_annotated_variants[4]
    assert len(v.parsed_csq) == 5
    for csq in v.parsed_csq:
        assert csq["GMAF"] == "T:0.0008"
        assert csq["ExAC_MAF"] == "T:0.0006"
        assert csq["ExAC_Adj_MAF"] == "T:8.649e-05"
