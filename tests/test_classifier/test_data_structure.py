from pathlib import Path

import pytest

from charger.classifier import CharGerResult
from charger.config import ACMG_MODULES, CHARGER_MODULES
from charger.variant import Variant


@pytest.fixture(scope="module")
def charger_result(test_root: Path):
    variant = next(
        Variant.read_vcf(
            test_root / "examples" / "grch37_vep85_5_variants.vcf", parse_csq=True
        )
    )
    return CharGerResult(variant)


def test_chargerresult_default_decision(charger_result: CharGerResult):
    acmg_modules = [m for _, modules in ACMG_MODULES.items() for m in modules]
    charger_modules = [m for _, modules in CHARGER_MODULES.items() for m in modules]

    assert set(charger_result.acmg_decisions) == set(acmg_modules)
    for module, decision in charger_result.acmg_decisions.items():
        assert decision is None

    assert set(charger_result.charger_decisions) == set(charger_modules)
    for module, decision in charger_result.charger_decisions.items():
        assert decision is None
