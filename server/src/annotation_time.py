import time
from config import DATA_DIR
from sys import stderr 

def start_time(username,docname):
	f = open(DATA_DIR+'/'+username+'/'+docname+'_time.txt','a+')
	f.write('START: ' + str(time.time())+'\n')
	f.close()
	return {}

def end_time(username,docname):
        f = open(DATA_DIR+'/'+username+'/'+docname+'_time.txt','a+')
        f.write('END: ' + str(time.time())+'\n')
        f.close()
	return {}
