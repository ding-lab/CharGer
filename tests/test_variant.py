import pytest

from charger.variant import AnnotatedVariant


@pytest.fixture
def snp():
    return AnnotatedVariant("1", 42, 42, "A", "G", {})


@pytest.fixture
def insertion():
    return AnnotatedVariant("1", 42, 42, "A", "ATAT", {})


@pytest.fixture
def deletion():
    return AnnotatedVariant("1", 42, 44, "ATA", "A", {})


@pytest.fixture
def sv_del_like():
    return AnnotatedVariant("1", 2827693, 2827701, "CCCCTCGCA", "C", {"SVTYPE": "DEL"})


def test_is_snp(snp, insertion, deletion, sv_del_like):
    assert snp.is_snp
    assert not insertion.is_snp
    assert not deletion.is_snp
    assert not sv_del_like.is_snp


def test_is_sv(snp, sv_del_like):
    assert not snp.is_sv
    assert sv_del_like.is_sv


def test_is_indel(snp, insertion, deletion, sv_del_like):
    assert not snp.is_indel
    assert insertion.is_indel
    assert deletion.is_indel
    assert not sv_del_like.is_indel
