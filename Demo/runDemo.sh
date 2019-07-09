bp1="BP1GeneList.txt"
pp2="PP2GeneList.txt"
out="charged.demo.tsv"
in="demo.vcf"
pvs1="inheritanceGeneList.txt"
pm5="somaticHotspots.hotspot3d.clusters"
#Test CharGer using BioMine ReST VEP & ClinVar
## grch37 VEP endpoint currently down as of Jan 2018
##charger -f demo.vcf -o charged.demo.tsv -D --inheritanceGeneList inheritanceGeneList.txt --PP2GeneList PP2GeneList.txt --BP1GeneList BP1GeneList.txt -l -E
#Test CharGer using BioMine ReST VEP & ClinVar
charger -f ${in} -o ${out} -D --inheritanceGeneList ${pvs1}
