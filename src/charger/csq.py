from collections import UserDict
from typing import Any, Dict, List, Set

from loguru import logger

logger.disable("charger")  # Disable emit logs by default

ALL_CONSEQUENCE_TYPES: List[str] = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]
"""All the possible consequence types fetched from `Ensembl v99`_ (January 2020).

The consequence types here are ordered by their severeness.

.. _Ensembl v99: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
"""

ALL_TRUNCATION_TYPES: List[str] = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
]
"""All consequence types considered as a truncation."""

ALL_INFRAME_TYPES: List[str] = [
    "inframe_insertion",
    "inframe_deletion",
    "stop_lost",
]
"""All consequence types considered as to be inframe."""


class CSQ(UserDict):
    """
    Consequence of a variant. Access each CSQ field like a `dict`.

    The class is used to set the annotation records in a :class:`~charger.variant.Variant`
    object. List of CSQ per feature will be stored at :attr:`Variant.parsed_csq
    <charger.variant.Variant.parsed_csq>`.

    Examples:

        >>> csq = variant.parsed_csq[0]; csq
        CSQ(SYMBOL='FANCM', HGVSc='ENST00000267430.5:c.5101N>T', Consequence='stop_gained', …)
        >>> list(csq.keys())[:5]
        ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene']
        >>> list(csq.values())[:5]
        ['T', 'stop_gained', 'HIGH', 'FANCM', 'ENSG00000187790']
        >>> csq['HGVSc']
        'ENST00000267430.5:c.5101N>T'
    """

    data: Dict[str, Any]

    __slots__ = [
        "data",
    ]

    #: Required CSQ fields. Will raise a `ValueError` if any of the fields is
    #: missing when creating a new CSQ object.
    REQUIRED_FIELDS: Set[str] = set(
        [
            "Allele",
            "Consequence",
            "SYMBOL",
            "Gene",
            "Feature_type",
            "Feature",
            "BIOTYPE",
            "HGVSc",
            "HGVSp",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "Existing_variation",
            "STRAND",
        ]
    )

    def __init__(self, dict=None, **kwargs):
        super().__init__(dict, **kwargs)
        missing_fields = self.REQUIRED_FIELDS - set(self.data.keys())
        if missing_fields:
            raise ValueError(
                f"CSQ misses these required fields: {', '.join(missing_fields)}"
            )

    @property
    def consequence_types(self) -> List[str]:
        """Get all the consequence types separated."""
        return self.data["Consequence"].split("&")

    def rank_consequence_type(self) -> int:
        """Rank the severeness of its consequence type (CSQ column ``Consequence``).

        Severe consequence type has smaller rank (smallest being 0). Ranking is based on the
        order in :attr:`ALL_CONSEQUENCE_TYPES`. When the CSQ has multiple consequence types
        separated by ``&``, return the smallest rank of all the types. When the consequence type
        is not known, return the biggest possible rank + 1.
        """
        ranks: List[int] = []
        for ct in self.consequence_types:
            try:
                rank = ALL_CONSEQUENCE_TYPES.index(ct)
            except ValueError:
                # Assign unknown consequence type to the lowest rank
                rank = len(ALL_CONSEQUENCE_TYPES)
                logger.warning(
                    "Got unknown consequence type: {ct}; assign its rank = {rank}",
                    ct=ct,
                    rank=rank,
                )
            ranks.append(rank)
        return min(ranks)

    def is_truncation_type(self) -> bool:
        """Whether the consequence type is truncation.

        See :attr:`ALL_TRUNCATION_TYPES` for the full list of consequence types.
        """
        return any(ct in ALL_TRUNCATION_TYPES for ct in self.consequence_types)

    def is_inframe_type(self) -> bool:
        """Whether the consequence type is inframe.

        See :attr:`ALL_INFRAME_TYPES` for the full list of consequence types.
        """
        return any(ct in ALL_INFRAME_TYPES for ct in self.consequence_types)

    def __repr__(self):
        fields = ["SYMBOL", "HGVSc", "Consequence"]
        details = []
        for f in fields:
            details.append(f"{f}={self.data[f]!r}")
        return f"CSQ({', '.join(details)}, …)"
