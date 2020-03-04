0.6.0 (*unreleased*)
====================
Backward-incompatible Changes
-----------------------------
- Rewrite CharGer to be Python 3.6+ compatible and Python 2.7 is no longer supported. Most of the options have been renamed and unified.
- Input variants must be normalized and biallelic, and must be annotated by VEP with the recommended parameters. Only VCF format is accepted (``.vcf``, ``.vcf.gz``, or ``.bcf``).
- ClinVar variant table must be Tabix indexed.
- v0.6.0 results are identical to v0.5.4 except for the conditions listed below:

  - ClinVar matching now requires the alternative allele to be identical.


0.5.4 (2019-09-30)
==================
Changes
-------
- Add support for vcf files annotated with VEP ≥ 90.
- When present in VEP annotation, gnomAD population frequencies are prioritized over ExAC.


0.5.3 (2019-09-30)
==================
Bug Fixes
---------
- Fix parsing ClinVar information.
- Fix ``parseMacPathogenicity()`` to handle variants with multiple submitters that receive both benign and pathogenic classifications, but no conflict is reported (i.e. ``isPathogenic == 1 and isBenign == 1 and isConflicted == 0``).
- [:pr:`19`] Fix ``var.splitHGVSc`` doesn't consider strand information. When ``override = True``, the ref and alt for genomic variants will be wrongly changed for minus strand transcripts. Default to ``override = False``.
- [:pr:`19`] Fix coordinates in ``getMacClinVarTSV()`` to match ``readVCF()``.
