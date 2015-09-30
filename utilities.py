def write_data_to_file( file, data, append=True ):
	if append:
		handle = open( file, 'a')
	else:
		handle = open( file, 'w')
	handle.write(data)
	handle.close()

def write_params_to_file( file, params, extra=None, header=True, append=True, format="%.3f" ):
	if append:
		handle = open( file, 'a')
	else:
		handle = open( file, 'w')
	keys = params.keys()
	if header:
		handle.write("\t".join(keys))
		handle.write("\n")
	values = [format%(params[k]) for k in keys]
	handle.write("\t".join(values))
	if extra:
		handle.write(extra)
	handle.write("\n")
	handle.close()

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
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')