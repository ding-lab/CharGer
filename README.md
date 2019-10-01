# CharGer

CharGer (Characterization of Germline variants) is a software tool for interpreting and predicting clinical pathogenicity of germline variants. CharGer gathers evidence from databases and annotations, provided by local tools and files or via ReST APIs, and classifies variants according to ACMG guidelines for assessing variant pathogenicity. User-designed pathogenicity criteria can be incorporated into CharGer’s flexible framework, thereby allowing users to create a customized classification protocol.

If you use CharGer, please cite our publication so we can continue to support CharGer development:

Adam D Scott, Kuan-Lin Huang, Amila Weerasinghe, R Jay Mashl, Qingsong Gao, Fernanda Martins Rodrigues, Matthew A Wyczalkowski, Li Ding, CharGer: clinical Characterization of Germline variants, Bioinformatics, Volume 35, Issue 5, 01 March 2019, Pages 865–867, https://doi.org/10.1093/bioinformatics/bty649

### Requirements
* python 2.7.x
* pip 10.x
* virtualenv (RECOMMENDED; assumed below)
* git / wget / unzip / curl (depending on the approach taken)

### Standard installation

(1) Set up a python virtual environment:
```sh
 mkdir -p /path/to/workdir
 cd /path/to/workdir
 virtualenv mycharger --python=python2.7
 cd mycharger
 . bin/activate
```

(2) Prepare for CharGer
```sh
 pip --version
```
If the indicated version of pip is < 10.x, you will first need to upgrade your pip because python.org has ended its support for the TLSv1.0 and TLSv1.1 protocols:
```sh
 curl https://bootstrap.pypa.io/get-pip.py | python
```

(3) Select one of the following installation methods:

* Binary modules option (i.e., the easy approach via <a href="https://pypi.org">PyPI</a>)
```sh
  pip install charger
```
&nbsp; &nbsp; &nbsp; &nbsp;This command downloads and installs CharGer and its dependencies. The charger executable is placed into your mycharger/bin directory and should be ready for use. Proceed to the Run section below.

* Source code option

&nbsp; &nbsp; &nbsp; &nbsp; Download the CharGer source using one of the following:
```sh
  git clone https://github.com/ding-lab/CharGer.git
```
&nbsp; &nbsp; &nbsp; &nbsp; or
```sh
  wget -O CharGer.zip https://github.com/ding-lab/CharGer/archive/master.zip
  unzip CharGer.zip
  mv CharGer-master CharGer
```
&nbsp; &nbsp; &nbsp; &nbsp;Then install CharGer and its dependencies:

```sh
  cd CharGer
  pip install .

  # Update your PATH environment variable
  # It is suggested also to append this line to your ~/.bash_profile or ~/.bashrc
  export PATH="/path/to/workdir/mycharger/CharGer/bin:${PATH}"
```

## Installation using conda

In case you do not have python 2.7.x installed in your machine or the installation process above does not run smoothly for you, please follow the installation steps below using conda:

(1) Download and install Anaconda® 

Download the installer for Anaconda® for python 2.7.x here (https://www.anaconda.com/download/) or download and install via command line as follows (please copy link for the .sh file appropriate for your operational system in the webpage provided):

```sh
  wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
```

Install anaconda:

```sh
  bash Anaconda2-5.2.0-Linux-x86_64.sh
```

(2) Prepare for CharGer

Create a virtual environment:
```sh
  conda create --name CharGer python=2.7
```

Switch to virtual environment:
```sh
  conda activate CharGer
```

Ensure that the appropriate version of pip (>10.x) is available in your environment:
```sh
  pip --version
```

If the indicated version of pip is < 10.x, you will first need to install or upgrade your pip because python.org has ended its support for the TLSv1.0 and TLSv1.1 protocols. 
To install pip in your environment, run:

```sh
  conda install pip
```

To upgrade pip in your environment, run:

```sh
  conda update pip
```

Install pysam in your environment:

```sh
  conda install pysam
```

(3) Download and install CharGer

```sh
  wget -O CharGer.zip https://github.com/ding-lab/CharGer/archive/master.zip
  unzip CharGer.zip
  mv CharGer-master/ CharGer
```

Install CharGer:

```sh
  cd CharGer
  pip install .
```

(4) Update your PATH environment variable

After installation is complete, leave the conda environment and update your PATH environment variable. 

Leave conda environment:

```sh
  conda deactivate
```
Update PATH (it is suggested that you also append this line to your ~/.bash_profile or ~/.bashrc) :
  
```sh
  export PATH=“path/to/anaconda2/envs/CharGer/bin:${PATH}”
```

## Run
Example for a VCF file

	charger -f <variant file> -o <output file>

## For Help
To obtain a summary of options and default values, type

	charger

For other help/support, please submit an issue with us <a href="https://github.com/ding-lab/CharGer/issues">here</a>.


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
	--include-vcf-details (flag)
Name your output file; otherwise it will be called charger_summary.tsv. 
You can opt to make the output into an HTML page, instead of a readable .tsv. 
If you need to be assured of properly linked URL's, use the url test flag. 

### Access data
	-l ClinVar (flag)
	-x ExAC (flag)
	-E VEP (flag)
	-t TCGA cancer types (flag)
Using these flags turns on accession features built in. 
For the ClinVar, ExAC, and VEP flags, if no local VEP or database is provided, then BioMine will be used to access the ReST interface. CharGer is currently capable of handling all VEP releases up until release 97. 
The TCGA flag allows disease determination from sample barcodes in a .maf when using a diseases file (see below). 

### Suppress data or overrides
	-O override with ClinVar description (flag)
	-D suppress needing disease specific (flag)
You can have CharGer override its pathogenic characterization with whatever ClinVar has. 
Suppressing disease specific variants takes any variants in the diseases file (see below) and treats them as equally pathogenic without disease consideration.

### Cross-reference data files
	-z pathogenic variants, .vcf
	-e expression matrix file, .tsv
	--inheritanceGeneList inheritance gene list file, (format: gene\tdisease\tmode_of_inheritance) .txt
	--PP2GeneList PP2 gene list file, (format: column of genes) .txt
	--BP1GeneList BP1 gene list file, (format: column of genes) .txt
	-d diseases file, (format: gene\\tdisease\\tmode_of_inheritance) .tsv
	-n de novo file, standard .maf
	-a assumed de novo file, standard .maf
	-c co-segregation file, standard .maf
	-H HotSpot3D clusters file, .clusters
Variants or genes from each of these files can be used as additional known information. 
An expression matrix file has columns for each sample, and its rows are genes. 
The genes should be approved HUGO symbols. 
HotSpot3D clusters can be used for versions v1.x.x. 

### Thresholds
  --recurrence-threshold HotSpot3D recurrence threshold (default = 2)
  --rare-threshold Allele frequency threshold for rare (default = 0.0005 (0.05%)):
  --common-threshold Allele frequency threshold for common (default = 0.005 (0.5%)):
The recurrence threshold will be pulled from the recurrence/weight column of the .clusters file when provided.

### Pathogenicity/benignity standard modules and scores
Specify the option and positive whole number value to change the default value.

Standard modules:
```
  --PVS1 very strong pathogenicity (default = 8)
  --PS1 , --PS2 , --PS3 , --PS4 strong pathogenicity (defaults: PS1 = 7, PS2=PS3=PS4 = 4)
  --PM1 , --PM2 , --PM3 , --PM4 , --PM5 , --PM6 moderate pathogenicity (defaults: PM1=PM2=PM3=PM4=PM5 = 2)
  --PP1 , --PP2 , --PP3 , --PP4 , --PP5 supporting pathogenicity (defaults: PP1=PP2=PP3=PP4=PP5 = 1)
  --BP1 , --BP2 , --BP3 , --BP4 , --BP5 , --BP6 , --BP7 supporting benignity (defaults: BP1=BP2=BP3=BP4=BP5=BP6=BP7 = -1)
  --BS1 , --BS2 , --BS3 , --BS4 strong benignity (defaults: BS1=BS2=BS3=BS4 = -4)
  --BA1 stand-alone benignity (defaults: BA1 = -8)
```

CharGer-defined modules and scores
```
  --PSC1 strong pathogenicity (defaults: PSC1 = 4)
  --PMC1 moderate pathogenicity (defaults: PMC1 = 2)
  --PPC1 , --PPC2 supporting pathogenicity (defaults: PPC1=PPC2 = 1)
  --BMC1 moderate benignity (defaults: BMC1 = -2)
  --BSC1 strong benignity (defaults: BSC1 = -6)
```

### Pathogenicity/benignity category thresholds
Specify the option and positive whole number value to change the default value.
```
  --min-pathogenic-score threshold for classifying variant as pathogenic (default = 9)
  --min-likely-pathogenic-score threshold for classifying variant as likely pathogenic (default = 5)
  --max-likely-benign-score threshold for classifying variant as likely benign (default = -4)
  --max-benign-score threshold for classifying variant as benign (default = -8)
```

### Local VEP
	--perl Path to Perl
	--vep-script Path to VEP
	--vep-config config-file for VEP
	--vep-cache Path to VEP cache directory
	--vep-version VEP version (default = 87)
	--vep-output VEP output file (default = charger.vep.vcf)
	--grch assembly GRCh verion (default = 37)
	--ensembl-release Ensembl release version (default = 75)
	--reference-fasta VEP reference fasta
	--fork Number of forked processes used in VEP (default = 0) 
This currently only works with .vcf input only. 
CharGer is currently capable of handling all VEP releases up until release 97. 
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

### Local databases (suppresses ReST)
	--exac-vcf ExAC vcf.gz
	--mac-clinvar-tsv ClinVar from MacArthur lab (clinvar_alleles.tsv.gz)
Using local databases suppresses the BioMine accession too. 
These files can be downloaded from their respective sites.

### Filters
	--frequency-filter Keep if allele frequency is lower (default = 1.0, process variant with any frequency):
	--vcf-any-filter Keep variants that do not pass all filters in .vcf input (flag)
	--mutation-types Keep types, as a comma-delimited list (no spaces)
Using filters will limit the variants processed. 
The rare option takes variants with allele frequency less than the given value. 
The .vcf any filter accepts only variants that have passed all filters. 
If no .vcf pass filter status given, the .vcf null value will be taken as having passed. 
Mutation types filtering requires a comma delimitted list (no spaces) using terms from Ensembl's consequence terms.

### ReST batch sizes
	-v VEP number of variants (default/max allowed = 300)
	-b ClinVar summary number of variants (default/max allowed = 500)
	-B ClinVar searchsize number of variants (default/max allowed = 50)
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
Use these for .tsv and/or .maf input variant files to specify columns of relevant data. 
CharGer makes use of genomic and protein variant annotations, so the more data made available, the better your results.
