from contextlib import closing
from pathlib import Path
from typing import Generator, Type, TypeVar

from cyvcf2 import VCF
from cyvcf2 import Variant as CyVCF2Variant

# A type hint variable to annotated the Factory method
# See https://github.com/python/typing/issues/58 for details
V = TypeVar("V", bound="AnnotatedVariant")


class AnnotatedVariant:
    """
    Representation of VEP annotated VCF biallelic variant.

    For normal usage, consider using :py:meth:`~read_vcf` to construct the objects from a VEP annotated VCF.

    Examples:
    """

    @classmethod
    def _from_cyvcf2(cls: Type[V], variant: CyVCF2Variant) -> V:
        """
        Create one AnnotatedVariant object from one
        :py:class:`cyvcf2.Variant <cyvcf2.cyvcf2.Variant>` VCF record.
        """
        return cls()

    @classmethod
    def read_vcf(cls: Type[V], path: Path) -> Generator[V, None, None]:
        """
        Create AnnotatedVariant object per VCF record from `path`.

        This function iterates the given VCF using :py:class:`cyvcf2.VCF <cyvcf2.cyvcf2.VCF>`.

        >>> vcf_reader = AnnnotatedVariant.read_vcf('my.vcf')
        >>> next(vcf_reader)
        """
        with closing(VCF(str(path))) as vcf:
            for variant in vcf:
                yield cls._from_cyvcf2(variant)
