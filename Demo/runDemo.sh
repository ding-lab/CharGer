#Test CharGer using BioMine ReST VEP & ClinVar
../bin/charger -f demo.vcf -o charged.demo.tsv -D --inheritanceGeneList inheritanceGeneList.txt --PP2GeneList PP2GeneList.txt --BP1GeneList BP1GeneList.txt -l -E
#Test CharGer using local VEP & ClinVar
#TODO
#../bin/charger -f demo.vcf -o charged.demo.tsv -D --inheritanceGeneList inheritanceGeneList.txt --PP2GeneList PP2GeneList.txt --BP1GeneList BP1GeneList.txt
