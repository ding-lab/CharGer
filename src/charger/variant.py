import re
from contextlib import closing
from enum import Enum, Flag, auto
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Type, TypeVar

import attr
from cyvcf2 import VCF
from cyvcf2 import Variant as CyVCF2Variant
from loguru import logger

from .csq import CSQ

logger.disable("charger")  # Disable emit logs by default


# A type hint variable to annotated the Factory method
# See https://github.com/python/typing/issues/58 for details
V = TypeVar("V", bound="Variant")


@attr.s(auto_attribs=True, repr=False, eq=False, order=False, slots=True)
class Variant:
    """
    Biallelic variant.

    For normal usage, consider using :meth:`~read_vcf` to construct the objects from a VEP annotated VCF.

    Examples:

        >>> variant = Variant('13', 32340300, 32340301, 'GT', 'G', id='rs80359550')
        >>> variant
        Variant(13:32340300GT>G info: )
        >>> v.is_snp()
        False
        >>> v.is_sv()
        False
        >>> v.is_indel()
        True
        >>> v.is_deletion()
        True

        Annotate it with online VEP,

            >>> v = next(Variant.read_vcf('rs80359550.vcf', parse_csq=True))
            >>> v
            Variant(13:32340300GT>G info: CSQ[4 parsed])
    """

    chrom: str  #: Chromosome.
    start_pos: int  #: Start position (1-based closed). Same as POS in the VCF record.
    end_pos: int  #: End position (1-based closed).
    ref_allele: str  #: Reference allele sequence.
    alt_allele: str  #: Alternative allele sequence (currently only allow one possible allele).
    id: Optional[str] = None
    """ID in the VCF record. `None` when the original value is ``.``."""
    filter: Optional[List[str]] = None
    """FILTER in the VCF record. `None` when the original value is ``PASS``."""
    info: Dict[str, Any] = attr.Factory(dict)
    """INFO in the VCF record."""

    parsed_csq: Optional[List["CSQ"]] = None
    """
    All parsed CSQ annotations of the variant as a list of :class:`.CSQ` objects.
    Use :meth:`read_vcf(parse_csq=True) <read_vcf>` to automatically parse CSQ
    while reading an annotated VCF.
    """

    _most_servere_csq: Optional["CSQ"] = attr.ib(default=None, init=False)
    """Cached most servere CSQ based on the consequence type."""

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

    def is_snp(self) -> bool:
        """`True` if the variant is a SNP."""
        if len(self.ref_allele) > 1:
            return False
        elif self.alt_allele not in ("A", "C", "G", "T"):
            return False
        return True

    def is_sv(self) -> bool:
        """`True` if the variant ia an SV."""
        return "SVTYPE" in self.info and self.info["SVTYPE"] is not None

    def is_indel(self) -> bool:
        """`True` if the variant ia an INDEL."""
        is_sv = self.is_sv()
        if len(self.ref_allele) > 1 and not is_sv:
            return True

        if self.alt_allele == ".":
            return False
        elif len(self.alt_allele) != len(self.ref_allele) and not is_sv:
            return True
        return False

    def is_deletion(self) -> bool:
        """`True` if the variant is a deletion."""
        if not self.is_indel():
            return False
        else:
            return len(self.ref_allele) > len(self.alt_allele)

    def get_most_servere_csq(self) -> "CSQ":
        """Get the most severe CSQ based on the consequence type.

        If multiple CSQs have the same consequence type, the canonical CSQ determined by VEP will be selected.
        """
        # Return the cache the result if available
        if self._most_servere_csq is not None:
            return self._most_servere_csq

        if self.parsed_csq is None:
            raise ValueError(
                "Variant {self!r} may not have annotation or its CSQ is not parsed."
            )

        rank_and_canonical_per_csq = []
        for csq_ix, csq in enumerate(self.parsed_csq):
            # Severe consequence have smaller rank
            rank_consequence_type = csq.rank_consequence_type()
            # Canonical transcript wins whenever there is a tie in the ranks of consequence type.
            try:
                if csq["CANONICAL"] == "YES":
                    rank_canonical = 0
                else:
                    rank_canonical = 1
            except KeyError:
                rank_canonical = 1
            rank_and_canonical_per_csq.append(
                (rank_consequence_type, rank_canonical, csq_ix)
            )
        rank_ct, rank_canonical, csq_ix = min(rank_and_canonical_per_csq)
        self._most_servere_csq = self.parsed_csq[csq_ix]
        return self._most_servere_csq

    def _parse_csq(self, csq_fields: List[str]):
        """Parse the CSQ info string based on the CSQ field spec.

        It returns a list of consequences per annotation(transcript)
        as a list of :class:`.CSQ` objects.
        """
        all_csq = self.info["CSQ"].split(",")
        parsed_csq_per_annotation = []
        for one_csq in all_csq:
            csq_values = one_csq.split("|")
            if len(csq_values) != len(csq_fields):
                raise ValueError(
                    f"Number of CSQ values (n={len(csq_values)}) and columns (n={len(csq_fields)})"
                    f"don't match. columns: {csq_fields!r}, values: {csq_values!r}"
                )
            parsed_csq_per_annotation.append(CSQ(dict(zip(csq_fields, csq_values))))
        self.parsed_csq = parsed_csq_per_annotation
        self.info["CSQ"] = self.parsed_csq

    def __repr__(self):
        type_name = type(self).__name__
        ref = limit_seq_display(self.ref_allele)
        alt = limit_seq_display(self.alt_allele)

        info_keys = set(self.info.keys())
        if self.parsed_csq is not None:
            info_keys.discard("CSQ")
            info_keys.add(f"CSQ[{len(self.parsed_csq)} parsed]")
        info_repr = f"info: {','.join(info_keys)}"

        return f"{type_name}({self.chrom}:{self.start_pos}{ref}>{alt} {info_repr})"

    @classmethod
    def from_cyvcf2(cls: Type[V], variant: CyVCF2Variant) -> V:
        """
        Create one Variant object based on the given
        :class:`cyvcf2.Variant <cyvcf2.cyvcf2.Variant>` VCF record.
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

    @classmethod
    def get_vep_csq_fields(cls: Type[V], vcf: VCF):
        """Extract the CSQ fields VEP output in the given VCF."""
        # Reverse the header order because the newer header appears later
        vcf_raw_headers = list(reversed(vcf.raw_header.splitlines()))
        # Find VEP version
        try:
            vep_header = next(l for l in vcf_raw_headers if l.startswith("##VEP="))
            vep_version = re.match(r"^##VEP=['\"]?v(\d+)['\"]?", vep_header).group(1)  # type: ignore
        except (StopIteration, AttributeError):
            logger.warning(f"Cannot find VEP version in the VCF header")
            vep_version = "UNKNOWN"

        # Get CSQ spec
        try:
            csq_info_header = next(
                l for l in vcf_raw_headers if l.startswith("##INFO=<ID=CSQ,")
            )
        except StopIteration:
            raise ValueError(f"Cannot find CSQ format in the VCF header")
        m = re.search(r"Format: ([\w\|]+)['\"]", csq_info_header)
        if m:
            csq_format = m.group(1)
        else:
            raise ValueError(
                f"Cannot parse the CSQ field format from its INFO VCF header: {csq_info_header}"
            )
        csq_fields = csq_format.split("|")

        logger.debug(
            f"VEP version {vep_version} with CSQ format [{len(csq_fields)} fields]: {','.join(csq_fields)}"
        )
        return csq_fields

    @classmethod
    def read_vcf(
        cls: Type[V], path: Path, parse_csq: bool = False
    ) -> Generator[V, None, None]:
        """
        Create one object per VCF record from `path`.

        This function walks through each variant record in the given VCF
        using :class:`cyvcf2.VCF <cyvcf2.cyvcf2.VCF>`,
        and yields the record as a :class:`Variant` object.

        Args:
            path: Path to the VCF.
            parse_csq: whether to parse the VEP annotated CSQ annotations.
                If `True`, the parsed CSQ will be stored in the generated
                :attr:`Variant.parsed_csq`.

        Returns:
            An generator walking through all variants per record.

        Examples:

            Read an annotated VCF::

                >>> vcf_reader = Variant.read_vcf('my.vcf', parsed_csq=True)
                >>> variant = next(vcf_reader)
                >>> variant
                Variant(14:45658326C>T info: CSQ[5 parsed])
                >>> variants[4].parsed_csq[0]
                CSQ(SYMBOL='FANCM', HGVSc='ENST00000267430.5:c.5101N>T', Consequence='stop_gained', …)

            Iterate all the VCF variants records::

                >>> for variant in vcf_reader:
                ...     print(variant.chrom, variant.parsed_csq[0]['Allele'])
        """
        with closing(VCF(str(path))) as vcf:
            if parse_csq:
                # Get the CSQ field definition from the VCF header
                csq_fields = cls.get_vep_csq_fields(vcf)
            for cy_variant in vcf:
                variant = cls.from_cyvcf2(cy_variant)
                if parse_csq:
                    variant._parse_csq(csq_fields)
                yield variant


def limit_seq_display(seq: str, limit: int = 5) -> str:
    """Limit the display of the sequence.

    Examples:

    >>> limit_seq_display('ATATCCG')
    'ATATC…'
    >>> limit_seq_display('ATA')
    'ATA'
    >>> limit_seq_display('ATA', limit=1)
    'A…'
    """
    if len(seq) > limit:
        seq = seq[:limit] + "…"
    return seq


class GeneInheritanceMode(Flag):
    """All possible modes of the gene inheritance dominance.

    Used by :attr:`CharGerConfig.inheritance_gene_table
    <charger.config.CharGerConfig.inheritance_gene_table>`.
    """

    AUTO_DOMINANT = auto()  #: The gene is autosomal dominant.
    AUTO_RECESSIVE = auto()  #: The gene is autosomal recessive.
    X_LINKED_DOMINANT = auto()  #: The gene is X-linked dominant.
    X_LINKED_RECESSIVE = auto()  #: The gene is X-linked recessive.
    Y_LINKED = auto()  #: The gene is Y-linked.

    @classmethod
    def parse(
        cls: Type["GeneInheritanceMode"], value: str
    ) -> Optional["GeneInheritanceMode"]:
        """Parse the inheritance modes from the given string. Multiple modes are comma separated.

        >>> m = GeneInheritanceMode.parse("autosomal dominant, autosomal recessive")
        >>> m
        <GeneInheritanceMode.AUTO_RECESSIVE|AUTO_DOMINANT: 3>

        >>> bool(m & GeneInheritanceMode.AUTO_RECESSIVE)
        True
        >>> bool(m & GeneInheritanceMode.Y_LINKED)
        False

        >>> GeneInheritanceMode.parse("unknown") is None
        True
        """
        MODE_TO_FLAG = {
            "autosomal dominant": cls.AUTO_DOMINANT,
            "autosomal recessive": cls.AUTO_RECESSIVE,
            "x-linked dominant": cls.X_LINKED_DOMINANT,
            "x-linked recessive": cls.X_LINKED_RECESSIVE,
            "y-linked": cls.Y_LINKED,
        }
        mode_flags = []
        # Convert the input value to be lower case only
        for mode in value.lower().split(","):
            mode = mode.strip()
            if mode == "unknown":
                # Unknown has no flag
                continue
            try:
                mode_flags.append(MODE_TO_FLAG[mode])
            except KeyError as e:
                raise ValueError(
                    f"Invalid variant inheritance mode {mode} from {value}"
                ) from e
        if mode_flags:
            # Combined all the flags
            combined_mode = mode_flags[0]
            for flag in mode_flags[1:]:
                combined_mode |= flag
            return combined_mode
        else:
            return None


class ClinicalSignificance(Enum):
    """All possible clinical significance types of a variant."""

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"
    UNCERTAIN = "Uncertain Significance"

    @classmethod
    def parse_clinvar_record(
        cls: Type["ClinicalSignificance"], record: Dict[str, str]
    ) -> "ClinicalSignificance":
        """Determine the pathogenicity of a ClinVar record."""
        # The default pathogenicity
        clin_sig = cls.UNCERTAIN

        is_pathogenic = (
            int(record["pathogenic"]) > 0 or int(record["likely_pathogenic"]) > 0
        )
        is_benign = int(record["benign"]) > 0 or int(record["likely_benign"]) > 0

        if int(record["conflicted"]):
            # Conflicted record has uncertain clinical significance
            return cls.UNCERTAIN

        # There may be multiple assertions
        all_clin_sig_assertions = record["clinical_significance"].lower().split("/")
        if is_benign and is_pathogenic:
            # Fix parsing of conflicting ClinVar assertion
            # by checking the ClinVar assertion
            for assertion in all_clin_sig_assertions:
                if "likely benign" in assertion:
                    clin_sig = cls.LIKELY_BENIGN
                elif "likely pathogenic" in assertion:
                    clin_sig = cls.LIKELY_PATHOGENIC
                elif "benign" in assertion:
                    clin_sig = cls.BENIGN
                    break
                elif "pathogenic" in assertion:
                    clin_sig = cls.PATHOGENIC
                    break
        elif is_benign:
            for assertion in all_clin_sig_assertions:
                if "likely benign" in assertion:
                    clin_sig = cls.LIKELY_BENIGN
                elif "benign" in assertion:
                    clin_sig = cls.BENIGN
                    break
        elif is_pathogenic:
            for assertion in all_clin_sig_assertions:
                if "likely pathogenic" in assertion:
                    clin_sig = cls.LIKELY_PATHOGENIC
                elif "pathogenic" in assertion:
                    clin_sig = cls.PATHOGENIC
                    break
        return clin_sig
