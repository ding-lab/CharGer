0.6.0 (*unreleased*)
====================
Improvements
------------


Bug Fixes
---------


0.5.4 (2019-09-30)
==================
Improvements
------------
- Added support for vcf files annotated with VEP releases ≥ 90.
- When present in VEP annotation, CharGer prioritizes gnomAD population frequencies over ExAC.


0.5.3 (2019-09-30)
==================
Improvements
------------

Bug Fixes
---------
- Fixed parsing ClinVar information.
- Fixed ``parseMacPathogenicity()`` to handle variants with multiple submitters that received both benign and pathogenic classifications, but no conflict is reported (i.e. ``isPathogenic == 1 and isBenign == 1 and isConflicted == 0``).
- [:pr:`19`] Fixed ``var.splitHGVSc`` doesn't consider strand information. When ``override = True``, the ref and alt for genomic variants would be wrongly changed for minus strand transcripts. Changed to ``override = False``.
- [:pr:`19`] Fixed coordinates in ``getMacClinVarTSV()`` to match ``readVCF()``.
