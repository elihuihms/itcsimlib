"""Utilities for file input/output.


"""

import os
import pickle
import scipy
import scipy.special
import scipy.signal

from collections import OrderedDict

from . import __version__


def write_data_to_file( file, data, append=True ):
	if append:
		handle = open( file, 'a')
	else:
		handle = open( file, 'w')
	handle.write(data)
	handle.close()

def write_params_to_file( file, params, append=True, header=True, pre=None, post=None, format="%.5E"):
	if append:
		handle = open( file, 'a')
	else:
		handle = open( file, 'w')

	if header:
		handle.write("\t".join( params.keys()))
		handle.write("\n")
		
	if pre:
		handle.write(str(pre))
		handle.write("\t")
	
	values = [format%(params[k]) for k in params.keys()]
	handle.write("\t".join(values))
	
	if post:
		handle.write("\t")
		handle.write(str(post))
		
	handle.write("\n")
	handle.close()
	
def read_params_from_file( file, row=1, header=0 ):
	counter,Hline,Dline = 0,None,None
	with open( file, 'r' ) as f:
		for l in f:
			if counter == row:
				Dline = l
			if counter == header:
				Hline = l
			counter+=1
	if Hline == None:
		raise Exception("Specified header row # (%i) not found in %s"%(row,file))
	if Dline == None:
		raise Exception("Specified data row # (%i) not found in %s"%(row,file))
	head,data = Hline.split(),map(float,Dline.split())
	return OrderedDict(zip(head,data))
	
def read_itcsimlib_exp( file, exp_args={} ):
	from .itc_experiment import ITCExperiment

	ignore = ("itcsim","Date","Ivol","units")
	data,h = scipy.genfromtxt(file,unpack=True),open(file)
	kwargs = {'Cell':{},'Syringe':{}}
	for a in [l.split()[1:] for l in h.readlines() if l[0]=='#']:
		if a == [] or a[0] in ignore:
			continue
		elif a[0] == 'Cell' or a[0] == 'Syringe':
			kwargs[a[0]][a[1]] = float(a[2])
		elif a[0].lower() == 'skip':
			kwargs['skip'] = list(map(int,a[1:]))
		else:
			kwargs[a[0]] = float(a[1])
	h.close()
	if not 'title' in kwargs:
		kwargs['title'] = os.path.splitext(os.path.basename(file))[0]
	
	# overwrite any file-obtained info with explicit values
	kwargs.update(exp_args)
	if len(data) == 2:
		return ITCExperiment(injections=data[0],dQ=data[1],**kwargs)
	elif len(data) == 3:
		return ITCExperiment(injections=data[0],dQ=data[1],dQ_err=data[2],**kwargs)
	else:
		return None # TODO : parser errors

def read_itcsimlib_pickle( path ):
	with open(path, 'rb') as handle:
		experiment =  pickle.load(handle)
		del experiment._itcsimlib_version # TODO : version checking
		return experiment

def write_itcsimlib_pickle( path, experiment ):
	experiment._itcsimlib_version = __version__
	with open(path, 'wb') as handle:
		return pickle.dump(experiment, handle, protocol=pickle.HIGHEST_PROTOCOL)

def read_nitpic_exp( file, exp_args={}, recalc_concs=False ):
	from .itc_experiment import ITCExperimentBase

	with open( file, 'rb' ) as buf:
		nitpic = pickle.load(buf, fix_imports=True, encoding="latin1", errors="strict")

	# TODO: more sanity checks
	assert len(nitpic['inj_vols']) == len(nitpic['NDH'])
	
	# get uncertainty in dh by scaling NDH interval (high and low, assuming +/- 1 SD) against NDH
	DH, NDH = scipy.array(nitpic['dh']), scipy.array(nitpic['NDH'])
	NDH_u, NDH_d = scipy.array(nitpic['NDHerrorsUp']), scipy.array(nitpic['NDHerrorsDown'])

	kwargs = {
		'skip' : [],
		'dQ_err' : DH * ( (NDH_u+NDH_d)/(2.0*NDH) ),
		'cellRef' : "M", 'syringeRef' : "X",
		'units' : 'cal' ,
		'title' : nitpic['inputFilename']}

	kwargs.update(exp_args)
	experiment = ITCExperimentBase(
		nitpic['experimental_temp'],
		nitpic['cell_V'],
		nitpic['inj_vols'],
		nitpic['dh'],
		{'M':nitpic['CellConc']},
		{'X':nitpic['SyrConc']},
		**kwargs)

	if not recalc_concs: # do we overwrite the concentrations we calculated with those in the file?
		for i in range(experiment.npoints):
			experiment.Concentrations[i]['M'] = nitpic['Mt'][i]
			experiment.Concentrations[i]['X'] = nitpic['Xt'][i]

	experiment.initialized = True

	return experiment
	
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
	return scipy.signal.savgol_filter(y, window_length=window_size, polyorder=order, deriv=deriv, delta=rate)	
