from typing import Dict, List, Tuple

import pytest
from pysam import TabixFile

from charger.classifier import CharGer
from charger.variant import Variant


@pytest.fixture(scope="module")
def clinvar_tabix(test_root):
    return TabixFile(
        str(
            test_root.joinpath("examples/annotations/clinvar_chrom_22_only.b37.tsv.gz")
        ),
        encoding="utf8",
    )


@pytest.fixture(scope="module")
def matched_variants() -> List[Tuple[Variant, Dict[str, str]]]:
    return [
        # snp
        (
            Variant("22", 29083907, 29083907, "G", "A", id="VCV000422506.2"),
            {"symbol": "CHEK2", "hgvs_c": "NM_007194.4:c.1610C>T"},
        ),
        # deletion
        (
            Variant("22", 50967717, 50967719, "AGG", "A", id="VCV000223079.1"),
            {"symbol": "TYMP", "hgvs_c": "NM_001953.4:c.263_264del"},
        ),
        # indel
        (
            Variant("22", 51065766, 51065767, "GA", "AG", id="VCV000003062.1"),
            {"symbol": "ARSA", "hgvs_c": "NM_001085428.2:c.34_35delTCinsCT"},
        ),
    ]


@pytest.fixture(scope="module")
def not_found_variants() -> List[Variant]:
    return [
        Variant("22", 29083907, 29083907, "G", "T"),
        Variant("22", 50967718, 50967719, "GG", "G"),
        Variant("22", 51065766, 51065767, "GA", "C"),
    ]


def test_match(clinvar_tabix, matched_variants: List[Tuple[Variant, Dict[str, str]]]):
    cols = clinvar_tabix.header[0][len("#") :].split("\t")
    for variant, known_clinvar_record in matched_variants:
        record = CharGer._match_clinvar_one_variant(variant, clinvar_tabix, cols)
        assert record is not None
        for known_col, known_val in known_clinvar_record.items():
            assert record[known_col] == known_val


def test_not_found(clinvar_tabix, not_found_variants: List[Variant]):
    cols = clinvar_tabix.header[0][len("#") :].split("\t")
    for variant in not_found_variants:
        record = CharGer._match_clinvar_one_variant(variant, clinvar_tabix, cols)
        assert record is None


def test_warn_ref_mismatch(clinvar_tabix, caplog):
    # wrong ref should still match the clinvar but with a warning
    cols = clinvar_tabix.header[0][len("#") :].split("\t")
    v_wrong_ref = Variant("22", 29083907, 29083907, "C", "A")
    record = CharGer._match_clinvar_one_variant(v_wrong_ref, clinvar_tabix, cols)
    assert record is not None
    assert (
        "got a clinvar match but their reference alleles are different" in caplog.text
    )
