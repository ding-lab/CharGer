#!/bin/python
#03 February 2016 - Kuan-Lin Huang @ WashU - 
 
import sys
import getopt

def main():
    def usage():
        print """
    liftover_CharGer_result.py : why do I exist?

    USAGE: liftover_CharGer_result.py [-h] <filename>
     -h    print this message
     <filename>    input file
        """

    #use getopt to get inputs
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h') #:after option meaning required arguments
    except getopt.GetoptError:
        print "liftover_CharGer_result.py <charGerFile> <inputFile>"

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
    
    CharGerHeader = inFile.readline()
    varCharGer = {}
    #read input file
    for line in charGerF:
        line=line.strip()
        F = line.split("\t")
        var = F[0]
        varCharGer[var]=line
    charGerF.close()


    try:
        inputF = open(args[1],"r")
    except IOError:
        print("File , args[1], does not exist!")

    header = inFile.readline()
    print header

    #read input file
    for line in inputF:
        line=line.strip()
        F = line.split("\t")
        chrom = F[0]
        start = F[1]
        stop = F[2]
        ref = F[3]
        alt = F[4]
        sample = F[21]
        genomeVar = str(sample) + "::" \
        + str(chrom) + ":" \
        + str(start) + "-" \
        + str(stop) \
        + str(ref) \
        + ">" + str(alt)
        
        CharGerAnno = "NA\tNA\tNA\tNA"
        if genomeVar in varCharGer:
            CharGerAnno = varCharGer[genomeVar]
        
        print "\t".join(F[0:42]) + "\t" + CharGerAnno  
        #print line + "\t" + CharGerAnno
    
    inputF.close()

if __name__ == "__main__":
    main()
