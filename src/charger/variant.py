from contextlib import closing
from pathlib import Path
from typing import Any, Dict, Generator, Optional, Type, TypeVar

from cyvcf2 import VCF
from cyvcf2 import Variant as CyVCF2Variant

# A type hint variable to annotated the Factory method
# See https://github.com/python/typing/issues/58 for details
V = TypeVar("V", bound="AnnotatedVariant")


class AnnotatedVariant:
    """
    Representation of VEP annotated VCF biallelic variant.

    For normal usage, consider using :py:meth:`~read_vcf` to construct the objects from a VEP annotated VCF.

    Args:
        chrom: Chromosome
        start_pos: Start position (1-based closed)
        end_pos: End position (1-based closed)
        ref_allele: Reference allele sequence
        alt_allele: Alternative allele sequence (the variant must be biallelic)

    Examples:
    """

    def __init__(
        self, chrom, start_pos, end_pos, ref_allele, alt_allele, id=None, raw_info=None
    ):
        self.chrom: str = chrom
        self.start_pos: int = start_pos
        self.end_pos: int = end_pos
        self.ref_allele: str = ref_allele
        self.alt_allele: str = alt_allele
        self.id: Optional[str] = id
        if raw_info is None:
            raw_info = {}
        self.raw_info: Dict[str, Any] = raw_info

        # Variant must be normalized
        if ref_allele == ".":
            raise ValueError(
                "ref_allele cannot be missing ('.'). Try normalize the variant first."
            )
        if alt_allele is None or alt_allele == ".":
            raise ValueError(
                "alt_allele cannot be missing ('.' or None). Try normalize the variant first."
            )

    __slots__ = [
        "chrom",
        "start_pos",
        "end_pos",
        "ref_allele",
        "alt_allele",
        "id",
        "raw_info",
    ]

    @property
    def is_snp(self) -> bool:
        """True if the variant is a SNP."""
        if len(self.ref_allele) > 1:
            return False
        elif self.alt_allele not in "ACGT":
            return False
        return True

    @property
    def is_sv(self) -> bool:
        """True if the variant ia an SV."""
        return "SVTYPE" in self.raw_info and self.raw_info["SVTYPE"] is not None

    @property
    def is_indel(self) -> bool:
        """True if the variant ia an INDEL."""
        is_sv = self.is_sv
        if len(self.ref_allele) > 1 and not is_sv:
            return True

        if self.alt_allele == ".":
            return False
        elif len(self.alt_allele) != len(self.ref_allele) and not is_sv:
            return True
        return False

    @property
    def is_deletion(self) -> bool:
        """True if the variant is a deletion."""
        if not self.is_indel:
            return False
        else:
            return len(self.ref_allele) > len(self.alt_allele)

    @classmethod
    def _from_cyvcf2(cls: Type[V], variant: CyVCF2Variant) -> V:
        """
        Create one AnnotatedVariant object from one
        :py:class:`cyvcf2.Variant <cyvcf2.cyvcf2.Variant>` VCF record.
        """
        return cls(
            chrom=variant.CHROM,
            start_pos=variant.start + 1,
            end_pos=variant.end,
            ref_allele=variant.REF,
            alt_allele=variant.ALT[0],
            id=variant.id,
            raw_info=dict(variant.INFO),
        )

    @classmethod
    def read_vcf(cls: Type[V], path: Path) -> Generator[V, None, None]:
        """
        Create AnnotatedVariant object per VCF record from `path`.

        This function iterates the given VCF using :py:class:`cyvcf2.VCF <cyvcf2.cyvcf2.VCF>`.

        >>> vcf_reader = AnnnotatedVariant.read_vcf('my.vcf')
        >>> next(vcf_reader)
        """
        with closing(VCF(str(path))) as vcf:
            lineno = vcf.raw_header.count("\n")
            for lineno, variant in enumerate(vcf, start=lineno + 1):
                yield cls._from_cyvcf2(variant)
