# CharGer
Characterization of Germline variants
## Install

	pip install .

## Configure
Add the following to your PATH environment (add it in ~/.bash_profile or ~/.bashrc)

	export PATH="/path/to/charger/bin:${PATH}"

## Run
Example for a VCF file

	charger -f <variant file> -o <output file>

## For Help
Run:

	charger

For detailed help/support, email Adam:

	adamscott@wustl.edu

## Usage Details
### Input data
	-m Standard .maf
	-f Standard .vcf
	-T Custom .tsv
Variant data may be input via at least one variant file. 
This means that if variants are spread across several files, then you can input one of each type. 
For the .maf and .tsv, use the custom columns to determine which columns to use. 
Note that a standard .maf does not include protein annotations. 
Use the custom column for the peptide change column. 
If your .vcf has VEP annotations, then CharGer should be able to parse the information. 
This information will be added to your variants when available.

### Output
	-o output file
	-w output as HTML (flag)
	-k annotate input (flag)
	--run-url-test test url when creating links
Name your output file; otherwise it will be called charger_summary.tsv. 
You can opt to make the output into an HTML page, instead of a readable .tsv. 
If you need to be assured of properly linked URL's, use the url test flag. 

### Access data
	-l ClinVar (flag)
	-x ExAC (flag)
	-E VEP (flag)
	-t TCGA cancer types (flag)
Using these flags turns on accession features built in. 
For the ClinVar, ExAC, and VEP flags, if no local VEP or databse is provided, then BioMine will be used to access the ReST interface. 
The TCGA flag allows disease determination from sample barcodes in a .maf when using a diseases file (see below). 

### Suppress data or overrides
	-O override with ClinVar description (flag)
	-D suppress needing disease specific (flag)
You can have CharGer override its pathogenic characterization with whatever ClinVar has. 
Suppressing disease specific variants takes any variants in the diseases file (see below) and treats them as equally pathogenic without disease consideration.

### Cross-reference data
	-z pathogenic variants, .vcf
	-e expression matrix file, .tsv
	-g gene list file, (format: gene\\tdisease\\tmode_of_inheritance) .txt
	-d diseases file, (format: gene\\tdisease\\tmode_of_inheritance) .tsv
	-n de novo file, standard .maf
	-a assumed de novo file, standard .maf
	-c co-segregation file, standard .maf
	-H HotSpot3D clusters file, .clusters
	-r recurrence threshold (default = 2)
Variants or genes from each of these files can be used as additional known information. 
An expression matrix file has columns for each sample, and its rows are genes. 
The genes should be approved HUGO symbols. 
HotSpot3D clusters can be used for versions v1.x.x. 
The recurrence threshold will be pulled from the recurrence/weight column of the .clusters file when provided.

### Pathogenicity module scoring
Specify option and positive whole number value to change the default value.

Standard modules:
```
--PVS1 very strong pathogenicity (default = 1)
--PS1, --PS2, --PS3, --PS4 strong pathogenicity (defaults = 1)
--PM1, --PM2, --PM3, --PM4, --PM5, --PM6 moderate pathogenicity (defaults = 1)
--PP1, --PP2, --PP3, --PP4, --PP5 supporting pathogenicity (defaults = 1)
--BA1 stand-alone benignity (default = 1)
--BS1, --BS2, --BS3, --BS4 strong benignity (defaults = 1)
--BP1, --BP2, --BP3, --BP4, --BP5, --BP6, --BP7 supporting benignity (defaults = 1)
```

CharGer-defined modules:
```
--PSC1 strong pathogenicity (default = 1)
--PMC1 moderate pathogenicity (default = 1)
--PPC1, --PPC2 supporting pathogenicity (defaults = 1)
--BSC1 strong benignity (default = 1)
--BMC1 moderate benignity (default = 1)
```

### Pathogenicity category thresholds
Specify option and positive whole number value to change the default value.
```
--min-pathogenic-score threshold for classifying variant as pathogenic (default = 8)
--min-likely-pathogenic-score threshold for classifying variant as likely pathogenic (default = 5)
--min-benign-score threshold for classifying variant as benign (default = 8)
--min-likely-benign-score threshold for classifying variant as likely benign (default = 4)
```

### Local VEP
	--vep-script Path to VEP
	--vep-dir Path to VEP directory
	--vep-cache Path to VEP cache directory
	--vep-version VEP version (default = 87)
	--vep-output VEP output file (default = charger.vep.vcf)
	--grch assembly GRCh verion (default = 37)
	--ensembl-release Ensembl release version (default = 75)
	--reference-fasta VEP reference fasta
	--fork Number of forked processes used in VEP (default = 0) 
This currently only works with .vcf input only. 
Annotations are run with the VEP everything flag, so any local plugins will be used. 
The BioMine accession is also suppressed when using a local VEP installaltion. 
The VEP directory is not the same as would be given to VEP's --dir option. 
Instead it is the path to the directory with the VEP .pl file. 
The VEP script is the .pl file only. 
If not given, it will be /vep-dir/variant\_effect\_predictor.pl. 
The VEP cache directory is the same as would be given to VEP's --dir-cache option. 
If you have multiple VEP versions, then specify the version you want to use. 
This can be different from the Ensembl release option. 
VEP output is the same os would be given to VEP's -o option and should end with .vcf. 
The default output file will be called charger.vep.vcf. 
The GRCh reference genome can be set to either 37 or 38. 
The reference Fasta file will be deteremined automatically if not specified. 
If the reference Fasta file is constructed automatically, then if, for example, the VEP chache is ~/.vep/, the Ensembl release is 74, and the reference assembly is 37, then the reference Fasta file will be ~/.vep/homo\_sapiens/74\_GRCH37/Homo\_sapiens.GRCh37.74.dna.primary\_assembly.fa.gz.  

### Local databases
	--exac-vcf ExAC vcf.gz
	--mac-clinvar-tsv ClinVar from MacArthur lab (clinvar_alleles.tsv.gz)
Using local databases suppresses the BioMine accession too. 
These files can be downloaded from their respective sites.

### Filters
	--rare Allele frequency threshold for rare/common (default = 1, process variant with any frequency):
	--vcf-any-filter Allow variants that do not pass all filters in .vcf input (flag)
	--mutation-types Comma delimited list of types to allow
Using filters will limit the variants processed. 
The rare option takes variants with allele frequency less than the given value. 
The .vcf any filter accepts only variants that have passed all filters. 
If no .vcf pass filter status given, the .vcf null value will be taken as having passed. 
Mutation types filtering requires a comma delimitted list (no spaces) using terms from Ensembl's consequence terms.

### ReST batch sizes
	-v VEP (#variants, default/max allowed = 150)
	-b ClinVar summary (#variants, default/max allowed = 500)
	-B ClinVar searchsize (#variants, default/max allowed = 50)
ReST API's usually have limits on the amount of data sent or received. 
Exceeding these batch sizes would normally lead to warnings and/or IP blockage, but CharGer and BioMine try to keep batches at safe sizes. Last updated limits February 2017.

### Custom columns (0-based)
	-G HUGO gene symbol
	-X chromosome
	-S start position
	-P stop position
	-R reference allele
	-A alternate allele
	-s strand
	-M sample name
	-C codon
	-p peptide change
	-L variant classification
	-F allele frequency
Use these for .tsv and/or .maf input variant files to specify columns of relevant data. 
CharGer makes use of genomic and protein variant annotations, so the more data made available the better your results.
