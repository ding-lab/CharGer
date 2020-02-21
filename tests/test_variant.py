from pathlib import Path
from typing import List

import pytest

from charger.variant import (
    ClinicalSignificance,
    GeneInheritanceMode,
    Variant,
    limit_seq_display,
)


@pytest.fixture()
def snp() -> Variant:
    return Variant("1", 42, 42, "A", "G")


@pytest.fixture
def insertion() -> Variant:
    return Variant("1", 42, 42, "A", "ATAT")


@pytest.fixture
def deletion() -> Variant:
    return Variant("1", 42, 44, "ATA", "A")


@pytest.fixture
def sv() -> Variant:
    return Variant("1", 2827693, 2827701, "CCCCTCGCA", "C", info={"SVTYPE": "DEL"})


def test_variant_is_snp(
    snp: Variant, insertion: Variant, deletion: Variant, sv: Variant
):
    assert snp.is_snp()
    assert not insertion.is_snp()
    assert not deletion.is_snp()
    assert not sv.is_snp()


def test_variant_is_sv(
    snp: Variant, insertion: Variant, deletion: Variant, sv: Variant
):
    assert not snp.is_sv()
    assert not insertion.is_sv()
    assert not deletion.is_sv()
    assert sv.is_sv()


def test_variant_is_indel(
    snp: Variant, insertion: Variant, deletion: Variant, sv: Variant
):
    assert not snp.is_indel()
    assert insertion.is_indel()
    assert deletion.is_indel()
    assert not sv.is_indel()


def test_variant_is_deletion(
    snp: Variant, insertion: Variant, deletion: Variant, sv: Variant
):
    assert not snp.is_deletion()
    assert not insertion.is_deletion()
    assert deletion.is_deletion()
    assert not sv.is_deletion()


@pytest.fixture(scope="module")
def grch37_vep85_annotated_variants(test_root: Path):
    return list(
        Variant.read_vcf(
            test_root / "examples" / "grch37_vep85_5_variants.vcf", parse_csq=True
        )
    )


@pytest.fixture(scope="module")
def grch38_vep95_annotated_variants(test_root: Path) -> List[Variant]:
    return list(
        Variant.read_vcf(
            test_root.joinpath("examples/grch38_vep95_50_variants.vcf.gz"),
            parse_csq=True,
        )
    )


@pytest.fixture(scope="module")
def grch38_vep95_annotated_variants_info_fixed(test_root: Path) -> List[Variant]:
    return list(
        Variant.read_vcf(
            test_root.joinpath("examples/grch38_vep95_50_variants.info_fixed.vcf.gz"),
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
def test_read_vcf_detect_vep_version(
    test_root: Path, caplog, vcf_pth: Path, vep_ver: int
):
    reader = Variant.read_vcf(test_root.joinpath(vcf_pth), parse_csq=True)
    next(reader)
    assert f"VEP version {vep_ver}" in caplog.text


def test_read_vcf_grch37(grch37_vep85_annotated_variants: List[Variant]):
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


def test_read_vcf_grch38(grch38_vep95_annotated_variants: List[Variant]):
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


def test_variant_parse_csq_grch38(grch38_vep95_annotated_variants: List[Variant]):
    v = grch38_vep95_annotated_variants[42]
    assert v.parsed_csq is not None
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


def test_variant_parse_csq_grch37(grch37_vep85_annotated_variants: List[Variant]):
    parsed_csq = grch37_vep85_annotated_variants[0].parsed_csq
    assert parsed_csq is not None
    assert len(parsed_csq) == 11
    assert parsed_csq[10]["SYMBOL"] == "KLC3"
    assert parsed_csq[10]["ENSP"] == "ENSP00000466974"

    v = grch37_vep85_annotated_variants[4]
    assert v.parsed_csq is not None
    assert len(v.parsed_csq) == 5
    for csq in v.parsed_csq:
        assert csq["GMAF"] == "T:0.0008"
        assert csq["ExAC_MAF"] == "T:0.0006"
        assert csq["ExAC_Adj_MAF"] == "T:8.649e-05"


def test_limit_seq_display() -> None:
    assert limit_seq_display("ATATCCG") == "ATATC…"
    assert limit_seq_display("ATA") == "ATA"
    assert limit_seq_display("ATA", limit=1) == "A…"


def test_geneinheritancemode_parse():
    assert (
        GeneInheritanceMode.parse("X-linked recessive")
        == GeneInheritanceMode.X_LINKED_RECESSIVE
    )

    m = GeneInheritanceMode.parse("autosomal recessive, autosomal dominant, unknown")
    assert m & (GeneInheritanceMode.AUTO_DOMINANT | GeneInheritanceMode.AUTO_RECESSIVE)
    assert not m & GeneInheritanceMode.X_LINKED_DOMINANT


def test_geneinheritancemode_parse_unkonwn():
    assert GeneInheritanceMode.parse("unknown") is None
    assert GeneInheritanceMode.parse("unknown,Unknown") is None


def test_geneinheritancemode_parse_invalid():
    with pytest.raises(ValueError, match="Invalid variant inheritance mode"):
        GeneInheritanceMode.parse("y-linked recessive")

    with pytest.raises(ValueError, match="Invalid variant inheritance mode"):
        GeneInheritanceMode.parse("autosomal recessive, y-linked recessive")


@pytest.mark.parametrize(
    "clinvar_record,correct_result",
    [
        (
            ("0", "0", "2", "0", "0", "Likely benign"),
            ClinicalSignificance.LIKELY_BENIGN,
        ),
        (("0", "0", "0", "1", "0", "Benign"), ClinicalSignificance.BENIGN),
        (
            ("0", "0", "3", "2", "0", "Benign/Likely benign"),
            ClinicalSignificance.BENIGN,
        ),
        (("3", "0", "0", "0", "0", "Pathogenic"), ClinicalSignificance.PATHOGENIC),
        (
            ("0", "2", "0", "0", "0", "Likely pathogenic"),
            ClinicalSignificance.LIKELY_PATHOGENIC,
        ),
        (
            ("1", "1", "0", "0", "0", "Pathogenic/Likely pathogenic"),
            ClinicalSignificance.PATHOGENIC,
        ),
        (
            ("1", "0", "1", "0", "1", "Conflicting interpretations of pathogenicity"),
            ClinicalSignificance.UNCERTAIN,
        ),
        (
            ("0", "1", "1", "1", "1", "Conflicting interpretations of pathogenicity"),
            ClinicalSignificance.UNCERTAIN,
        ),
    ],
)
def test_clinical_significance_from_clinvar(clinvar_record, correct_result):
    record_cols = [
        "pathogenic",
        "likely_pathogenic",
        "likely_benign",
        "benign",
        "conflicted",
        "clinical_significance",
    ]
    record = dict(zip(record_cols, clinvar_record))
    assert ClinicalSignificance.parse_clinvar_record(record) is correct_result


@pytest.mark.parametrize(
    "clinvar_record,correct_result",
    [
        (
            ("1", "2", "1", "1", "0", "Likely benign"),
            ClinicalSignificance.LIKELY_BENIGN,
        ),
        (
            ("1", "2", "1", "1", "0", "Benign/Likely benign"),
            ClinicalSignificance.BENIGN,
        ),
        (
            ("1", "2", "1", "1", "0", "Likely pathogenic"),
            ClinicalSignificance.LIKELY_PATHOGENIC,
        ),
        (
            ("1", "2", "1", "1", "0", "Pathogenic/Likely pathogenic"),
            ClinicalSignificance.PATHOGENIC,
        ),
    ],
)
def test_clinical_significance_from_conflicting_clinvar_record(
    clinvar_record, correct_result
):
    record_cols = [
        "pathogenic",
        "likely_pathogenic",
        "likely_benign",
        "benign",
        "conflicted",
        "clinical_significance",
    ]
    record = dict(zip(record_cols, clinvar_record))
    assert ClinicalSignificance.parse_clinvar_record(record) is correct_result
