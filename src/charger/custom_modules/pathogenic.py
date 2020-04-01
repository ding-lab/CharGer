from typing import TYPE_CHECKING

from ..result import ModuleDecision
from ..variant import GeneInheritanceMode

if TYPE_CHECKING:
    # Only to import the type checking related modules during type check
    # to prevent circular import
    from ..classifier import CharGerResult, InheritanceGenesType


# region: Strong
def run_psc1(
    result: "CharGerResult", inheritance_genes: "InheritanceGenesType"
) -> None:
    """Run CharGer PSC1 module per variant."""
    most_severe_csq = result.variant.get_most_severe_csq()
    gene_symbol = most_severe_csq["SYMBOL"]
    if most_severe_csq.is_truncation_type() and gene_symbol in inheritance_genes:
        mode = inheritance_genes[gene_symbol]
        # Gene is autosomal dominant
        if mode is not None and mode & GeneInheritanceMode.AUTO_RECESSIVE:
            result.charger_decisions["PSC1"] = ModuleDecision.PASSED
            return
    result.charger_decisions["PSC1"] = ModuleDecision.FAILED


# endregion

# region: Moderate
def run_pmc1(
    result: "CharGerResult", inheritance_genes: "InheritanceGenesType"
) -> None:
    """Run CharGer PMC1 module per variant."""
    most_severe_csq = result.variant.get_most_severe_csq()
    if most_severe_csq.is_truncation_type():
        result.charger_decisions["PMC1"] = ModuleDecision.PASSED
    else:
        result.charger_decisions["PMC1"] = ModuleDecision.FAILED


# endregion

# region: Supporting
def run_ppc1(
    result: "CharGerResult", inheritance_genes: "InheritanceGenesType"
) -> None:
    """Run CharGer PPC1 module per variant."""
    most_severe_csq = result.variant.get_most_severe_csq()
    gene_symbol = most_severe_csq["SYMBOL"]
    if most_severe_csq.is_inframe_type() and gene_symbol in inheritance_genes:
        mode = inheritance_genes[gene_symbol]
        # Gene is autosomal dominant
        if mode is not None and mode & GeneInheritanceMode.AUTO_RECESSIVE:
            result.charger_decisions["PPC1"] = ModuleDecision.PASSED
            return
    result.charger_decisions["PPC1"] = ModuleDecision.FAILED


def run_ppc2(
    result: "CharGerResult", inheritance_genes: "InheritanceGenesType"
) -> None:
    """Run CharGer PPC2 module per variant."""
    most_severe_csq = result.variant.get_most_severe_csq()
    if most_severe_csq.is_inframe_type():
        result.charger_decisions["PPC2"] = ModuleDecision.PASSED
    else:
        result.charger_decisions["PPC2"] = ModuleDecision.FAILED


# endregion
