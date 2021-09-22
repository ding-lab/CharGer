from typing import TYPE_CHECKING, List

from ..result import ModuleDecision
from ..variant import GeneInheritanceMode, Variant

if TYPE_CHECKING:
    # Only to import the type checking related modules during type check
    # to prevent circular import
    from ..classifier import CharGerResult, InheritanceGenesType


# region: (Very) Strong
def run_pvs1(
    result: "CharGerResult", inheritance_genes: "InheritanceGenesType"
) -> None:
    """Run ACMG :term:`PVS1` module per variant."""
    # Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single
    # or multi-exon deletion) in a gene where LOF is a known mechanism of disease.
    #
    # Caveats:
    # - Beware of genes where LOF is not a known disease mechanism (e.g., GFAP, MYH7)
    # - Use caution interpreting LOF variants at the extreme 3′ end of a gene
    # - Use caution with splice variants that are predicted to lead to exon skipping but leave the
    #   remainder of the protein intact
    # - Use caution in the presence of multiple transcripts
    most_severe_csq = result.variant.get_most_severe_csq()
    gene_symbol = most_severe_csq["SYMBOL"]
    if most_severe_csq.is_truncation_type() and gene_symbol in inheritance_genes:
        mode = inheritance_genes[gene_symbol]
        # Gene is autosomal dominant
        if mode is not None and mode & GeneInheritanceMode.AUTO_DOMINANT:
            result.acmg_decisions["PVS1"] = ModuleDecision.PASSED
            # TODO: Check the expression effect if it's given
            return
    result.acmg_decisions["PVS1"] = ModuleDecision.FAILED


def run_ps1(result: "CharGerResult", pathogenic_variants: List[Variant]) -> None:
    """Run ACMG :term:`PS1` module per variant."""
    # Same amino acid change as a previously established pathogenic variant regardless of
    # nucleotide change.
    #
    # Caveats:
    # - Beware of changes that impact splicing rather than at the amino acid/protein level
    if result.clinvar is not None:
        pass


# endregion

# region: Moderate
def run_pm4(result: "CharGerResult", inheritance_genes: "InheritanceGenesType") -> None:
    """Run ACMG :term:`PM4` module per variant."""
    most_severe_csq = result.variant.get_most_severe_csq()
    gene_symbol = most_severe_csq["SYMBOL"]
    if most_severe_csq.is_inframe_type() and gene_symbol in inheritance_genes:
        mode = inheritance_genes[gene_symbol]
        # Gene is autosomal dominant
        if mode is not None and mode & GeneInheritanceMode.AUTO_DOMINANT:
            result.acmg_decisions["PM4"] = ModuleDecision.PASSED
            return
    result.acmg_decisions["PM4"] = ModuleDecision.FAILED


# endregion
