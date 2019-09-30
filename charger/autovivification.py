#!/usr/bin/env python
# autovivification - extends dict
# CharGer - Characterization of Germline variants
# author: 
#	- Adam D Scott (ascott@genome.wustl.edu)
#	- Fernanda Martins Rodrigues (fernanda@wustl.edu)
#	- Jay R. Mashl (rmashl@wustl.edu)
#	- Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.5.3 - September, 2019

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
