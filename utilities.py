import os
import os.path
import thermo

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

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
	
def read_itcsimlib_exp( file, **kwargs ):
	from scipy import genfromtxt
	ignore = ("itcsim","Date","Ivol","units")
	data,h = genfromtxt(file,unpack=True),open(file)
	info = {'Cell':{},'Syringe':{}}
	for a in [l.split()[1:] for l in h.readlines() if l[0]=='#']:
		if a == [] or a[0] in ignore:
			continue
		elif a[0] == 'Cell' or a[0] == 'Syringe':
			info[a[0]][a[1]] = float(a[2])
		elif a[0].lower() == 'skip':
			info['skip'] = map(int,a[1:])
		else:
			info[a[0]] = float(a[1])
	h.close()
	if not 'title' in info:
		info['title'] = os.path.splitext(os.path.basename(file))[0]
	return info,data
	
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
	# from http://stackoverflow.com/questions/22988882/how-to-smooth-a-curve-in-python
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) / 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')