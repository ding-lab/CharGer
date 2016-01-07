# CharGer
Characterization of Germline variants
## Dependencies
### Latest WebAPIs
CharGer is still under construction, so a one step install is not yet 
available. CharGer requires modules developed for several web API's developed 
by AdamDS. These can be installed by getting the latest commit of 
AdamDS/WebAPIs, found here
#### git clone https://github.com/AdamDS/WebAPIs.git
It may be necessary to periodically update the version of these modules, so if 
you already have the WebAPIs repository, in your local repository of WebAPIs, 
use the command:
#### git pull origin master
Once you have the latest WebAPIs, in the head directory of the repository, 
type:
#### python setup.py install
## Run
CharGer will (at some point) be able to run with only variant input OR with 
variant input and gene lists.
### Without gene lists
#### python charger.py -i <variant file> -o <output file>
