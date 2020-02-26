CharGer Internals
=================
Here documents CharGer's internal design by modules:

.. toctree::

   config
   classifier
   variant
   csq
   io
   console


Usage
-----
The module can be accessed externally by ``charger.<module_name>``. For example,

.. code-block::

   from charger.config import CharGerConfig
   from charger.classifier import CharGer
   from charger.variant import Variant
