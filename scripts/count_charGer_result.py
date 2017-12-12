#!/usr/bin/env python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt

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

def main():
    def usage():
        print """
    count_charGer_result.py : why do I exist?

    USAGE: count_charGer_result.py [-h] <filename>
     -h    print this message
     <filename>    input file
        """

    #use getopt to get inputs
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h') #:after option meaning required arguments
    except getopt.GetoptError:
        print "count_charGer_result.py  <charGer annotated File>"

    for opt, arg in opts: #store the input options
        if opt == '-h': # h means user needs help
            usage(); sys.exit()

    if len(args) < 1:
        usage(); sys.exit("input file missing")

    #open input file
    try:
        charGerF = open(args[0],"r")
    except IOError:
        print("File , args[0], does not exist!")
    
    cancer2gene2charger = autovivification({})
    cancer2gene2clinvar = autovivification({})

    #read input file
    for line in charGerF:
        line=line.strip()
        F = line.split("\t")
        gene = F[6]
        cancer = F[22]
        charger = F[44]
        clinvar = F[45]
        if charger in cancer2gene2charger[cancer][gene]:
            cancer2gene2charger[cancer][gene][charger] += 1
        else:
            cancer2gene2charger[cancer][gene][charger] = 1

        if clinvar in cancer2gene2clinvar[cancer][gene]:
            cancer2gene2clinvar[cancer][gene][clinvar] += 1
        else:
            cancer2gene2clinvar[cancer][gene][clinvar]=1

    charGerF.close()

    for cancer in cancer2gene2charger:
        for gene in cancer2gene2charger[cancer]:
            for charger in cancer2gene2charger[cancer][gene]:
                count = cancer2gene2charger[cancer][gene][charger]
                print '\t'.join( [ "CharGer", cancer, gene , charger , str(count) ] )

    for cancer in cancer2gene2clinvar:
        for gene in cancer2gene2clinvar[cancer]:
            for clinvar in cancer2gene2clinvar[cancer][gene]:
                count = cancer2gene2clinvar[cancer][gene][clinvar]
                print '\t'.join( [ "ClinVar", cancer, gene , clinvar , str(count) ] )


if __name__ == "__main__":
    main()
