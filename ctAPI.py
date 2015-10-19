#!/usr/bin/python
# author: Adam D Scott
# first created: 2015*10*07
##### Subset terminology
#	results - search results
#	show
##### Query terminology
#	term - general search terms, multiple terms joined by "+"
#	recr - recruitment status
#		Recruiting, 
#	no_unk - filter out studies with unknown status
#		Y, N
#	rslt - studies with/without results
#		With, Without
#	type - type of study
#		Intr, 
#	cond - conditions considered
#	intr - interventions used
#	titles - study titles
#	outc - measured outcomes
#	spons - sponsors
#	lead - lead sponsor
#	id - study ID
#	state1 (1-3) - state of study
#		NA%3AUS%3AAL = Alabama, state two letters following second %3A
#	cntry1 (1-3) - country of study
#		NA%3AUS = United States, ES%3ACN = China, 
#	locn - location of study
#	gndr - recruit gender
#		Female, Male
#	age - recruit age
#		0 = child (birth-17), 1 = adult (18-65), 2 = senior (66+)
#	phase - trial phase
#		phase0 = 4, phase1 = 0, phase2 = 1, phase3 = 2, phase4 = 3
#	fund - primary funder
#		0 = NIH, 1 = other US agency, 2 = industry, 3 = other
#	safe - safety issue outcome measures
#		Y, N
#	rcv_s - study first received start of range
#		MM%2FDD%2FYYYY
#	rcv_e - study first received end of range
#		MM%2FDD%2FYYYY
#	lup_s - study last updated start of range
#		MM%2FDD%2FYYYY
#	lup_e - study last updated end of range
#		MM%2FDD%2FYYYY
#
#	NOTE: %2C = ,
#https://clinicaltrials.gov/ct2/results?term=BRAF+cancer&recr=Recruiting&no_unk=Y&rslt=With&type=Intr&cond=cancer&intr=drug&titles=BRAF&outc=oume&spons=spco&lead=sple&id=stid&state1=NA%3AUS%3AAL&cntry1=NA%3AUS&state2=&cntry2=&state3=&cntry3=&locn=birmingham&gndr=Female&age=1&phase=0&phase=3&fund=0&safe=Y&rcv_s=04%2F11%2F1985&rcv_e=04%2F07%2F2000&lup_s=03%2F08%2F2003&lup_e=04%2F24%2F2010
#		
#https://clinicaltrials.gov/ct2/results?
###
#	endpoint	https://clinicaltrials.gov/ct2/
#	subset		results
#	query		cond=cancer&intr=drug
#	url			endpoint + subset + query

from restAPI import restAPI

class ctAPI(restAPI):
	endpoint = "https://clinicaltrials.gov/ct2/"
	results = "results"
	show = "show"
	def __init__(self,**kwargs):
		subset = kwargs.get("subset",'')
		if not subset:
			super(ctAPI,self).__init__(ctAPI.endpoint,ctAPI.results)
		else:
			if (subset == ctAPI.results or subset == ctAPI.show)
				super(ctAPI,self).__init__(ctAPI.endpoint,subset)
			else:
				print "ADSERROR: bad subset. restAPI.subset initializing to results"
				super(ctAPI,self).__init__(ctAPI.endpoint,ctAPI.results)
		self.url = ctAPI.endpoint + self.subset
	def __repr__(self):
		rep = "ctAPI\n\tendpoint = " + ctAPI.endpoint
		rep += "\n\turl = " + self.url
		return rep
	
	def setSubset(self,subset):
		self.subset = subset
		self.action = ""
		self.url = ctAPI.endpoint + subset
	def resetURL(self):
		self.action = ""
		self.url = ctAPI.endpoint + self.subset

	def beginQuery(self):
		self.action = "?"
