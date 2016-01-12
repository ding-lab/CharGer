#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value