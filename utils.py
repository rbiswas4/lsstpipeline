#!/usr/bin/env python



def get_userinfo():

	import socket
	import os 
	import time



	hostname =  socket.gethostname()
	username =  os.getlogin()

	string =  "Hostname: " + hostname + "\n"
	string += "Username: " + username + "\n"
	string += "Time    : "     + time.asctime() +"\n" 

	return string
def str2bool(string):

        string = 	string.lower()
	if string == "true" or string == "yes" or string == "1":
		return True
	elif string == "false" or string == "no"  or string == "0": 
		return False
	else :
		raise ValueError("Could not Interpret the boolean for %s" %string)

		
def readjob(string, dict):

	key, val  = string.split("=")


	filename  = None
	if ":" in val:
		vallist = val.split(":")
		value = vallist[0].strip()
		filename  = vallist[1].strip()

	else:
		value = val.strip()
	key  = key.strip()

	dict[key] = [value, filename] 

	return 0
	
def processfile(jobfile) :

	"""
	reads the list of jobs from a file jobfile to perform along with the 
	methods and setup for each job. 

	args: 
		jobfile : str, mandatory
			name of file containing the list of jobs to do 

	"""

	jobdict = {}

	f = open(jobfile) 
	line = f.readline()


	while line !="":
		#strip comments
		actual = line.split("#")[0] 
		if len(actual) != 0 :
			if "=" in actual:
				#print "The string is ", actual
				readjob(actual, jobdict)

		line = f.readline()

	f.close()
	return jobdict


if __name__=="__main__":


	jobs  = processfile ("config/joblist")

	print jobs
