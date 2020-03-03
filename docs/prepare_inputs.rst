.. prep_input:

Preparing the Inputs
====================
CharGer relies on the annotations in the input VCF, and additional annotations for different classification rules/modules. This page helps the user to prepare all the inputs for CharGer.

.. currentmodule:: charger.config

.. |CLI| replace:: :ref:`command line option <cli>`



Input variants
--------------
.. note::
    The resulting file goes to ``--input`` |CLI| or :attr:`CharGerConfig.input`.


Normalize the variants
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

    bcftools norm \
        -f /path/to/genome_ref.fa
        --multiallelics - \
        --check-ref e \
        -Oz -o input_normalized.vcf.gz \
        input.vcf.gz


Annotate the variants using Ensembl VEP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

    vep --format vcf --vcf \
        --assembly GRCh38 \
        --everything --af_exac \
        --offline \
        --cache --dir_cache /path/to/vep_cache/ \
        --fasta /path/to/Homo_sapiens.GRCh38.dna.toplevel.fa \
        --input_file input_normalized.vcf.gz \
        --output_file input_normalized.vep.vcf



Known pathogenic variants
-------------------------
.. note::
    The resulting file goes to ``--pathogenic-variant`` |CLI| or :attr:`CharGerConfig.pathognic_variant`.



Inheritance gene table
----------------------
.. note::
    The resulting file goes to ``--inheritance-gene-table`` |CLI| or :attr:`CharGerConfig.inheritance_gene_table`.



HotSpot3D results
-----------------
.. note::
    The resulting file goes to ``--hotspot3d-cluster`` |CLI| or :attr:`CharGerConfig.hotspot3d_cluster`.
