#!/usr/bin/python
# author: Adam D Scott (amviot@gmail.com)
# first created: 2015*09*28

#import requests
#import json
#from requests.auth import HTTPDigestAuth
#http://docs.python-requests.org/en/latest/user/authentication/?highlight=authentication
#CGT API section 2.3,
#"The API uses a process called Digest Authentication as a security
# measure to verify user requests submitted to the database."

class restAPI(object):
	'''REST API class, has 
		endpoint = the endpoint
		subset = the subset realm of the RESTful service
		action = the query, file upload, etc.'''
	def __init__(self,endpoint,subset):
		self.endpoint = endpoint
		self.subset = subset
		self.action = ""

	def queryAll(self):
		return self.site + "?query=NOT%20asdf"

	def partial(self,term,value):
		self.action += term + ":" + str(value)
	def addPartial(self,condition,term,value):
		self.action += "%20" + condition + "%20" + term + ":" + str(value)

	def exact(self,term,value):
		self.action += term + "=" + str(value)
	def addExact(self,condition,term,value):
		self.action += "%20" + condition + "%20" + term + "=" + str(value)

	def multiOr(self,term,values):
		values = map(str,values)
		stripped = [ v.strip() for v in values ]
		self.action += term + ":OR(" + ','.join(stripped) + ")"
	def addMultiOr(self,condition,term,values):
		values = map(str,values)
		stripped = [ v.strip() for v in values ]
		self.action += "%20" + condition + "%20" + term + ":OR(" + ','.join(stripped) + ")"

	def buildURL(self):
		return self.endpoint + self.subset + self.action

	def buildURLJSON(self):
		return self.endpoint + self.subset + self.action + "&_type=json"
