"""Basic thermodynamic transforms.


"""

import math


_R = 8.3144621 # J/(K*mol)
_UNITS = ('J','kJ','cal','kcal')

def convert_to_J(units,value):
	"""Convert a value in the specified units to Joules.
	
	Arguments
	---------
	units : string
		The units that the value parameter is in.
	value : float
		The value to convert to Joules.
	
	Returns
	-------
	float
		The value in Joules.
	"""
	assert units in _UNITS

	if units == 'J':
		return value
	elif units == 'kJ':
		return value*1000.0
	elif units == 'cal':
		return J_from_cal(value)
	elif units == 'kcal':
		return J_from_cal(value)*1000.0

def convert_from_J(units,value):
	"""Convert a value in Joules to the desired units.
	
	Arguments
	---------
	units : string
		The units that the desired value is to be returned in.
	value : float
		The value (in Joules) to convert.
	
	Returns
	-------
	float
		The value in the specified units.
	"""
	assert units in _UNITS
	
	if units == 'J':
		return value
	elif units == 'kJ':
		return value/1000.0
	elif units == 'cal':
		return cal_from_J(value)
	elif units == 'kcal':
		return cal_from_J(value)/1000.0	

def dG_from_Kd( Kd, T ):
	"""Convert a diassociation constant into a free energy change (in Joules) at a specified temperature.
	
	Arguments
	---------
	Kd : float
		A disassociation constant.
	T : float
		The desired temperature.
	
	Returns
	-------
	float
		The free energy change associated with the binding event at temperature T.
	"""
	return _R*T*math.log(Kd)

def Kd_from_dG( dG, T ):
	"""Convert a free energy change (in Joules) at a specified temperature to a diassociation constant.
	
	Arguments
	---------
	dG : float
		The free energy change of binding.
	T : float
		The desired temperature.
	
	Returns
	-------
	float
		The disassociation constant associated with the free energy change at temperature T.
	"""
	return math.exp( dG/(_R*T) )

def dS_from_dGdH( dG, dH, T ):
	"""Return the change in entropy for a process at some temperature from a given free energy and enthalpy change.
	
	Arguments
	---------
	dG : float
		The free energy change.
	dH : float
		The enthalpy change.
	T : float
		The desired temperature.
	
	Returns
	-------
	float
		The change in entropy.
	"""
	return (dG-dH)/(-1.0*T)

def J_from_cal( cal ):
	"""Convert an energy in calories to Joules.
	
	Arguments
	---------
	cal : float
		The energy value (in calories) to convert.
			
	Returns
	-------
	float
		The energy value, in Joules.
	"""
	return cal / 0.239005736

def cal_from_J( joules ):
	"""Convert an energy in Joules to calories.
	
	Arguments
	---------
	cal : float
		The energy value (in Joules) to convert.
			
	Returns
	-------
	float
		The energy value, in calories.
	"""
	return joules * 0.239005736

def get_ratios( A, B ):
	"""Return a list of ratios between the points in two lists, for non-numpy lists of numbers. (Why?)
	
	Arguments
	---------
	A : list of floats
		The numerators.
	B : list of floats
		The denominators.
	
	Returns
	-------
	list of floats
		The A/B ratios.
	"""
	n = len(A)
	assert( n == len(B) )
	ret = [0.0]*n
	for i in range(n):
		ret[i] = A[i] / B[i]
	return ret

def get_scale( A, B ):
	"""Return the average scaling factor between two lists of values.
	
	Arguments
	---------
	A : list of floats
		Curve or data set A.
	B : list of floats
		Curve or data set B.
	
	Returns
	-------
	float
		The average scaling value between A versus B
	"""
	n = len(A)
	assert( n == len(B) )
	return ( sum(A)/n )/( sum(B)/n )

def normalize( A, B ):
	"""Superimpose (normalize) the set of points B upon A.
	
	Arguments
	---------
	A : list of floats
		Curve or data set A.
	B : list of floats
		Curve or data set B.
	
	Returns
	-------
	list of floats
		B, superimposed on A
	"""
	n = len(A)
	assert( n == len(B) )

	# calculate root mean square for each dataset
	A_avg, B_avg = sum(A)/n, sum(B)/n
	A_rms, B_rms = 0.0, 0.0
	for i in range(n):
		A_rms += (A[i] - A_avg)**2
		B_rms += (B[i] - B_avg)**2
	A_rms, B_rms = math.sqrt(A_rms/n), math.sqrt(B_rms/n)

	# translate B data to origin, normalize scale, then transform back using data to be normalized to
	B_norm = [0.0]*n
	for i in range(n):
		B_norm[i] = (A_rms*(B[i] - B_avg)/B_rms) + A_avg

	return B_norm

def dK_Gibbs_Helmholtz( T, T0, K0, dH0, dCp ):
	"""Return the temperature-corrected change association constant from a reference association constant and reference enthalpy.
	
	Arguments
	---------
	T : float
		The temperature at which to return the corrected free energy.
	T0 : float
		The reference temperature.
	K0 : float
		The reference association constant (association constant at temperature T0).
	dH0 : float
		The reference change in enthalpy (enthalpy change at temperature T0).
	dCp : float
		The change in heat capacity.
		
	Returns
	-------
	float
		The corrected association constant at temperature T.

	Notes
	-----
		Inverted Eqn. 19 from Winzor and Jackson (2006), also Naghibi et al., 1995
		Assumes a constant (temperature-independent) change in heat capacity (i.e. linear dH w.r.t. T)
	"""

	return math.exp(
		math.log(K0) + ( ((dH0-(T0*dCp)) / _R)*((1.0/T0)-(1.0/T)) ) + ((dCp/_R)*math.log(T/T0))
		)

def dH_vant_Hoff( dH0, dCp, T, T0 ):
	"""Return the temperature-corrected change in enthalpy from a reference enthalpy.
	
	Arguments
	---------
	dH0 : float
		The reference change in enthalpy (enthalpy change at temperature T0).
	dCp : float
		The change in heat capacity.
	T : float
		The temperature at which to return the corrected enthalpy.
	T0 : float
		The reference temperature.
		
	Returns
	-------
	float
		The corrected enthalpy at temperature T.
	
	Notes
	-----
		Integrated van't Hoff equation 12a from Prabhu & Sharp, AR Reviews (2005)
		Assumes a constant (temperature-independent) change in heat capacity (i.e. linear dH w.r.t. T)
	"""
	return dH0 + (dCp*(T-T0))


def dG_vant_Hoff( dG0, dH0, dCp, T, T0 ):
	"""Return the temperature-corrected change in free energy from a reference free energy and reference enthalpy.
	
	Arguments
	---------
	dG0 : float
		The reference change in free energy (free energy change at temperature T0)
	dH0 : float
		The reference change in enthalpy (enthalpy change at temperature T0).
	dCp : float
		The change in heat capacity.
	T : float
		The temperature at which to return the corrected free energy.
	T0 : float
		The reference temperature.
		
	Returns
	-------
	float
		The corrected free energy at temperature T.

	Notes
	-----
		Integrated van't Hoff equation 12c from Prabhu & Sharp, AR Reviews (2005)
		Assumes a constant (temperature-independent) change in heat capacity (i.e. linear dH w.r.t. T)
	"""
	dS0 = (dH0 - dG0) / T0
	return dH0 - (T*dS0) + (dCp*( (T-T0) - (T*math.log(T/T0)) ))

def dG_vant_Hoff_dH( dG0, dH, dCp, T, T0 ):
	"""
	Return the temperature-corrected change in free energy from a reference free energy and a given enthalpy.
	
	Arguments
	---------
	dG0 : float
		The reference change in free energy (free energy change at temperature T0)
	dH : float
		The change in enthalpy at temperature T (use dG_vant_Hoff() if using a reference enthalpy).
	dCp : float
		The change in heat capacity.
	T : float
		The temperature at which to return the corrected free energy.
	T0 : float
		The reference temperature.
		
	Returns
	-------
	float
		The corrected free energy at temperature T.
	
	Notes
	-----
		Integrate van't Hoff equation, but using dH instead of dH0
		Assumes a constant (temperature-independent) change in heat capacity (i.e. linear dH w.r.t. T)

		Derivation from eqn. 12c, the integrated van't Hoff equation, from Prabhu and Sharp (2005):
		dG = dH(T0) - T*dS(T0) + dCp[ (T-T0) - T*ln(T/T0) ]
		...
		Substitute definition for dS0 from Gibbs-Helmholtz eqn. (see Winzor and Jackson):
		dH(T0) = dG(T0) + T0*dS(T0)
		dS(T0) = ( dH(T0) - dG(T0) )/T0

		dG = dH(T0) - (T/T0)*[ dH(T0) - dG(T0) ] + dCp[ (T-T0) - T*ln(T/T0) ]
		...
		Substitute definition for dH0 from Prabhu and Sharp (2005);
		dH(T) = dH(T0) + dCp*(T-T0)
		dH(T0) = dH(T) - dCp*(T-T0)

		dG = dH(T) - dCp*(T-T0) - (T/T0)*[ dH(T) - dCp*(T-T0) - dG(T0) ] + dCp[ (T-T0) - T*ln(T/T0) ]
	"""

	return dH - dCp*(T-T0) - (T/T0)*( dH - dCp*(T-T0) - dG0 ) + dCp*( (T-T0) - T*math.log(T/T0) )

