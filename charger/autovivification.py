#!/usr/bin/env python
# autovivification - extends dict
# author: Kuan-lin Huang (khuang@genome.wustl.edu) & Adam D Scott (ascott@genome.wustl.edu)
# version: v0.0 - 2016*01*12

class autovivification(dict):
	'''Implementation of perl's autovivification feature.'''
	def __init__( self , *args , **kwargs ):
		super( autovivification , self ).__init__( *args , **kwargs )
		self.itemlist = super( autovivification , self ).keys()
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self, item)
		except KeyError:
			value = self[item] = type(self)()
			return value
