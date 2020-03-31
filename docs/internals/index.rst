CharGer Internals
=================
Here documents CharGer's internal design by modules:

.. toctree::

   config
   result
   classifier
   acmg_modules
   variant
   csq
   io
   console


Reuse CharGer modules
---------------------
The module can be accessed externally by ``charger.<module_name>``. For example,

.. code-block::

   from charger.classifier import CharGer
   from charger.config import CharGerConfig
   from charger.console import setup_logger
   from charger.variant import Variant

The following script programmatically calls CharGer:

.. literalinclude:: ../../scripts/debug_example.py
   :lines: 10-
