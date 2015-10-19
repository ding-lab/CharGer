import Pair
import Variant

mu1 = Variant.Variant(gene="BRAF",HGVSp="p.V600E",transcript="ENST0987")
mu2 = Variant.Variant(chromosome=7,start=12345,reference="A",variant="T")
p = Pair.Pair(mu1,mu2,0.2,"1UWH")

p.addDistance(0.22,"1UWJ")

p.addDistance(0.1,"BLAH")

p.printDistances("\t")

p.printPair('\t')

q = p.getShortest()

print(q[0] + ' ' + str( q[1] ) )
