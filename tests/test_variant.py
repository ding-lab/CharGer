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


def test_is_snp(snp, insertion, deletion, sv):
    assert snp.is_snp
    assert not insertion.is_snp
    assert not deletion.is_snp
    assert not sv.is_snp


def test_is_sv(snp, insertion, deletion, sv):
    assert not snp.is_sv
    assert not insertion.is_sv
    assert not deletion.is_sv
    assert sv.is_sv


def test_is_indel(snp, insertion, deletion, sv):
    assert not snp.is_indel
    assert insertion.is_indel
    assert deletion.is_indel
    assert not sv.is_indel


def test_is_deletion(snp, insertion, deletion, sv):
    assert not snp.is_deletion
    assert not insertion.is_deletion
    assert deletion.is_deletion
    assert not sv.is_deletion


def test_read_vcf_grch37():
    variants = list(
        Variant.read_vcf(
            HERE / "examples" / "grch37_vep85_5_variants.vcf", parse_csq=True
        )
    )
    assert len(variants) == 5
    assert variants[0].chrom == "19"
    assert variants[0].start_pos == 45855804
    assert variants[0].end_pos == 45855805
    assert variants[0].ref_allele == "CT"
    assert variants[0].alt_allele == "C"
    assert variants[0].id is None
    assert variants[0].filter is None
    assert variants[0].is_indel
    assert variants[0].is_deletion
    assert not variants[0].is_snp
    assert not variants[0].is_sv
