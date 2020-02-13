import re
from contextlib import closing
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Type, TypeVar

import attr
from cyvcf2 import VCF
from cyvcf2 import Variant as CyVCF2Variant
from loguru import logger

logger.disable("charger")  # Disable emit logs by default

# A type hint variable to annotated the Factory method
# See https://github.com/python/typing/issues/58 for details
V = TypeVar("V", bound="Variant")


@attr.s(auto_attribs=True, repr=True, eq=False, order=False, slots=True)
class Variant:
    """
    Biallelic variant.

    For normal usage, consider using :py:meth:`~read_vcf` to construct the objects from a VEP annotated VCF.

    Args:
        chrom: Chromosome
        start_pos: Start position (1-based closed)
        end_pos: End position (1-based closed)
        ref_allele: Reference allele sequence
        alt_allele: Alternative allele sequence (the variant must be biallelic)

    Examples:
    """

    chrom: str
    start_pos: int
    end_pos: int
    ref_allele: str
    alt_allele: str
    id: Optional[str] = None
    filter: Optional[List[str]] = None
    info: Dict[str, Any] = attr.Factory(dict)

    def __attrs_post_init__(self):
        # Variant must be normalized
        if self.ref_allele == ".":
            raise ValueError(
                "ref_allele cannot be missing ('.'). Try normalize the variant first."
            )
        if self.alt_allele is None or self.alt_allele == ".":
            raise ValueError(
                "alt_allele cannot be missing ('.' or None). Try normalize the variant first."
            )

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
        return "SVTYPE" in self.info and self.info["SVTYPE"] is not None

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
    def from_cyvcf2(cls: Type[V], variant: CyVCF2Variant) -> V:
        """
        Create one object based on
        :py:class:`cyvcf2.Variant <cyvcf2.cyvcf2.Variant>` VCF record.
        """
        return cls(
            chrom=variant.CHROM,
            start_pos=variant.start + 1,
            end_pos=variant.end,
            ref_allele=variant.REF,
            alt_allele=variant.ALT[0],
            id=variant.ID,
            filter=variant.FILTER,
            info=dict(variant.INFO),
        )

    def _parse_csq(self, csq_columns: List[str]):
        csq_values = self.info["CSQ"].split("|")
        self.info["CSQ"] = dict(zip(csq_columns, csq_values))

    @classmethod
    def get_vep_csq_columns(cls: Type[V], vcf: VCF):
        vcf_raw_headers = vcf.raw_header.splitlines()
        # Find VEP version
        try:
            vep_header = next(
                l for l in reversed(vcf_raw_headers) if l.startswith("##VEP=")
            )
            vep_version = re.match(r"^##VEP=['\"]?v(\d+)['\"]?", vep_header).groups(1)  # type: ignore
        except (StopIteration, AttributeError):
            raise ValueError(f"Cannot find VEP information")

        # Get CSQ spec
        try:
            csq_info = vcf.get_header_type("CSQ")
            csq_format = re.search(r"Format: ([\w\|]+)['\"]$", csq_info["Description"]).group(1)  # type: ignore
            csq_columns = csq_format.split("|")
        except KeyError:
            raise ValueError(f"Cannot find CSQ INFO line")

        logger.debug(f"VEP version {vep_version} with CSQ format: {csq_format}")
        return csq_columns

    @classmethod
    def read_vcf(cls: Type[V], path: Path, parse_csq=False) -> Generator[V, None, None]:
        """
        Create one object per VCF record from `path`.

        This function iterates the given VCF using :py:class:`cyvcf2.VCF <cyvcf2.cyvcf2.VCF>`.

        >>> vcf_reader = Variant.read_vcf('my.vcf')
        >>> next(vcf_reader)
        """
        with closing(VCF(str(path))) as vcf:
            if parse_csq:
                csq_columns = cls.get_vep_csq_columns(vcf)
            lineno = vcf.raw_header.count("\n")
            for lineno, cy_variant in enumerate(vcf, start=lineno + 1):
                variant = cls.from_cyvcf2(cy_variant)
                if parse_csq:
                    variant._parse_csq(csq_columns)
                yield variant
