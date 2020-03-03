.. _cli:

Command line interface
======================
.. currentmodule:: charger.config


All options
-----------
Each options here has an one-to-one mapping to an attribute of :class:`CharGerConfig` by replacing all the dashes with underscores. For example, ``--input`` maps to :attr:`CharGerConfig.input`, and ``--pathogenic-variant`` maps to :attr:`CharGerConfig.pathogenic_variant`.

.. argparse::
    :ref: charger.console.create_console_parser
    :prog: charger
    :nodescription:
    :noepilog:


    --override-acmg-score : @after
        It overrides the default score for each ACMG module defined by :attr:`_default_acmg_scores` and stores the final scores at :attr:`CharGerConfig.acmg_module_scores`.

    --override-charger-score : @after
        It overrides the default score for each CharGer module defined by :attr:`_default_charger_scores` and stores the final scores at :attr:`CharGerConfig.charger_module_scores`.