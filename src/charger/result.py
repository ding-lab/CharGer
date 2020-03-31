from enum import Enum, auto
from typing import Any, Callable, Dict, List, Optional, Type

import attr

from .config import ACMG_MODULES, CHARGER_MODULES
from .variant import Variant


class ModuleDecision(Enum):
    """The decision of a ACMG or CharGer module of one variant.

    Used by :attr:`CharGerResult.acmg_decisions` and :attr:`CharGerResult.charger_decisions`.

    Examples:

        Skip PVS1 classification:

        >>> result = CharGerResult(Variant('19', 45855804, 45855804, 'CT', 'C'))
        >>> result.acmg_decisions['PVS1'] = ModuleDecision.SKIPPED
    """

    PASSED = auto()
    """The variant matched the criteria of the module."""
    FAILED = auto()
    """The variant didn't match the criteria of the module."""

    @classmethod
    def _gen_decision_template(
        cls: Type["ModuleDecision"], available_modules: Dict[str, List[str]]
    ) -> Callable[[], Dict[str, Optional["ModuleDecision"]]]:
        """
        Generate the decision template for all the available modules
        during the initiation of :class:`CharGerResult`.
        """

        def gen_template():
            decisions = {}
            for module_type, modules in available_modules.items():
                for m in modules:
                    decisions[m] = None
            return decisions

        return gen_template


@attr.s(auto_attribs=True, eq=False, order=False, slots=True)
class CharGerResult:
    """CharGer classification result of one variant.

    Examples:

        >>> variant = Variant(19, 45855804, 45855804, 'CT', 'C')
        >>> result = CharGerResult(variant)
        >>> result.acmg_decisions['PVS1'] is None
        True
    """

    variant: Variant
    """The input variant of this result."""

    clinvar: Dict[str, Any] = attr.Factory(dict)
    """ClinVar annotation of this variant."""

    acmg_decisions: Dict[str, Optional[ModuleDecision]] = attr.Factory(
        ModuleDecision._gen_decision_template(ACMG_MODULES)
    )
    """
    The classification decision per ACMG module of this variant.

    `None` if the module is not run. See :class:`ModuleDecision` for the possible decisions.
    """

    charger_decisions: Dict[str, Optional[ModuleDecision]] = attr.Factory(
        ModuleDecision._gen_decision_template(CHARGER_MODULES)
    )
    """
    The classification decision of each CharGer module of this variant.

    Same usage as :attr:`acmg_decisions`.
    """
