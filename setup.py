#https://docs.python.org/2/distutils/examples.html
from distutils.core import setup
version = "0.5.2"
setup( \
	name = 'CharGer' , 
	version = version , 
	author = 'Adam D Scott, Kuan-lin Huang, Amila Weerasinghe' ,
	author_email = 'adamscott@wustl.edu, kuan-lin.huang@wustl.edu, amila@wustl.edu' ,
	maintainer = 'Adam D Scott' ,
	maintainer_email = 'adamscott@wustl.edu' ,
	url = 'https://github.com/ding-lab/CharGer' ,
	description = 'Characterization of Germline variants' ,
	long_description = 'Characterization of Germline variants \
		Following ACMG guidelines for genomic variant pathogenicity \
		classification, CharGer cross-checks user variants against \
		several public databases for information. CharGer then \
		provides a pathogenicity score according to ACMG or \
		CharGer\'s custom scale.' ,
	download_url = 'https://github.com/ding-lab/CharGer/archive/v' + \
		version + '.tar.gz' ,
	scripts = ['bin/charger'] ,
	keywords = [ "bioinformatics" , "genomics" , "pathogenic" , "benign" , \
				"ClinVar" , "ExAC" , "VEP" , "BioMine" , "germline" , \
				"variant" , "classifier" ] ,
	classifiers = [ \
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)" ,
		"Programming Language :: Python" , 
		"Programming Language :: Python :: 2.7" , 
		"Development Status :: 4 - Beta" , 
		"Intended Audience :: Developers" , 
		"Intended Audience :: Science/Research" , 
		"Topic :: Internet" , 
		"Topic :: Scientific/Engineering" , 
		"Topic :: Scientific/Engineering :: Bio-Informatics" , 
		"Topic :: Scientific/Engineering :: Chemistry" , 
	] ,
	license = 'GPL-3.0' ,
	package_dir = { 'charger' : 'charger' , 
	} ,
	packages = [ \
		'charger' ,
	] , #each of the directories with modules (aka packages)
	install_requires = [ \
		'BioMine' ,
		'AdvancedHTMLParser' , 
		'requests' ,
		'pysam' , 
		'PyVCF' ,
		'SciPy' ,
		#'transvar' ,
	] , #auto installs with pip install
	#dependency_links = ['https://github.com/zwdzwd/transvar/archive/v2.1.23.20160321.tar.gz']
)
