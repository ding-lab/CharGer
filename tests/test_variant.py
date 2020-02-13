import pytest

from charger.variant import AnnotatedVariant


@pytest.fixture()
def snp(request):
    return AnnotatedVariant("1", 42, 42, "A", "G")


@pytest.fixture
def insertion():
    return AnnotatedVariant("1", 42, 42, "A", "ATAT")


@pytest.fixture
def deletion():
    return AnnotatedVariant("1", 42, 44, "ATA", "A")


@pytest.fixture
def sv():
    return AnnotatedVariant(
        "1", 2827693, 2827701, "CCCCTCGCA", "C", raw_info={"SVTYPE": "DEL"}
    )


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
