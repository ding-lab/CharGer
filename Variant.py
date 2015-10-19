class Variant(object):
	def __init__(self,**kwargs):
		self.gene = kwargs.get('gene',None)
		self.chromosome = kwargs.get('chromosome',None)
		self.start = kwargs.get('start',None)
		self.stop = kwargs.get('stop',None)
		self.reference = kwargs.get('reference',None)
		self.variant = kwargs.get('variant',None)
		self.HGVSp = kwargs.get('HGVSp',None)
		self.transcript = kwargs.get('transcript',None)
	def printVariant(self,delim):
		if self.gene != None:
			print self.gene + delim ,
		if self.chromosome != None:
			print str(self.chromosome) + delim ,
		if self.start != None:
			print str(self.start) + delim ,
		if self.stop != None:
			print str(self.stop) + delim ,
		if self.reference != None:
			print self.reference + delim ,
		if self.variant != None:
			print self.variant + delim ,
		if self.HGVSp != None:
			print self.HGVSp + delim ,
		if self.transcript != None:
			print self.transcript + delim ,
		print ""
	def attr(self):
		attributes = []
		if self.gene != None:
			attributes.append(self.gene)
		if self.chromosome != None:
			attributes.append(self.chromosome)
		if self.start != None:
			attributes.append(self.start)
		if self.stop != None:
			attributes.append(self.stop)
		if self.reference != None:
			attributes.append(self.reference)
		if self.variant != None:
			attributes.append(self.variant)
		if self.HGVSp != None:
			attributes.append(self.HGVSp)
		if self.transcript != None:
			attributes.append(self.transcript)
		return attributes

'''
mu = Variant(gene="BRAF",chromosome=7,start=12345,stop=123456,reference="AT",variant="GC",HGVSp="p.A123R")
mu.printVariant('\t')
nu = Variant(gene="BRAF",HGVSp="p.A123R")
nu.printVariant('\t')
ou = Variant(gene="BRAF",chromosome=7,start=12345,stop=123456,reference="AT",variant="GC")
ou.printVariant('\t')
pu = Variant(chromosome=7,start=12345,stop=123456,reference="AT",variant="GC")
pu.printVariant('\t')
qu = Variant(chromosome=7,start=12345,reference="AT",variant="GC")
qu.printVariant('\t')

print qu.attr()
'''
