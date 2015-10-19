import Variant

class Pair(object):
	numPairs = 0;
	def __init__(self,variant1,variant2,distance,pdb):
		self.variant1 = variant1
		self.variant2 = variant2
		self.distInfo = { pdb : distance }
		self.shortPDB = pdb
		self.shortDistance = distance
		Pair.numPairs += 1
	def __del__(self):
		Pair.numPairs -= 1
	def addDistance(self,newDist,newPDB):
		self.distInfo[newPDB] = newDist
		if self.shortDistance > newDist:
			self.shortDistance = newDist
			self.shortPDB = newPDB
		Pair.numPairs += 1
	def printPair(self,delim):
		print( delim.join( ( ':'.join( str(x) for x in self.variant1.attr() ) , ':'.join( str(x) for x in self.variant2.attr() ) , self.shortPDB , str( self.shortDistance ) ) ) )
	def printDistances(self,delim):
		print( delim.join( ( ':'.join( str(x) for x in self.variant1.attr() ) , ':'.join( str(x) for x in self.variant2.attr() ) , self.shortPDB , str( self.shortDistance ) ) ) )
		for pdb in self.distInfo:
			print( '\t' + pdb + ' ' + str( self.distInfo[pdb] ) )
	def getShortest(self):
		return (self.shortPDB,self.shortDistance)
