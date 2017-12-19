import os
import sys
import pdb
import unittest
#from biomine.variant.mafvariant import mafvariant
#from biomine.variant.clinvarvariant import clinvarvariant
from charger.chargervariant import chargervariant

class testchargervariant( unittest.TestCase ):
	#print( "asdfas" )
	def test_nonzero( self ):
		v = chargervariant()
		if ( v ):
			return True
		else:
			return False
if __name__ == '__main__':
	unittest.main()
