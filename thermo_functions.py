from math 	import exp,log

def get_ratios( A, B ):
	"""
	Returns a list of ratios between the points in two lists
	"""
	n = len(A)
	assert( n == len(B) )
	ret = [0.0]*n
	for i in xrange(n):
		ret[i] = A[i] / B[i]
	return ret

def get_scale( A, B ):
	"""
	Returns the correct scaling factor between two lists of values
	"""
	n = len(A)
	assert( n == len(B) )
	return ( sum(A)/n )/( sum(B)/n )

def dK_Gibbs_Helmholtz( R, T, T0, K0, dH0, dCp ):
	"""
	Eqn. 19 from Winzor and Jackson (2006), also Naghibi et al., 1995
	"""

	return exp(
		log(K0) + ( ((dH0-(T0*dCp)) / R)*((1.0/T0)-(1.0/T)) ) + ((dCp/R)*log(T/T0))
		)

def dH_vant_Hoff( dH0, dCp, T, T0 ):
	"""
	Integrated van't Hoff equation 12a from Prabhu & Sharp, AR Reviews (2005)
	"""
	return dH0 + (dCp*(T-T0))


def dG_vant_Hoff( dG0, dH0, dCp, T, T0 ):
	"""
	Integrated van't Hoff equation 12c from Prabhu & Sharp, AR Reviews (2005)
	Assumes a constant (i.e. linear) heat capacity
	"""
	dS0 = (dH0 - dG0) / T0
	return dH0 - (T*dS0) + (dCp*( (T-T0) - (T*log(T/T0)) ))

def dQ_calc( Q, V0, I_vol ):
	"""
	Eqn. 10 from Microcal's "ITC Data Analysis in Origin"
	"""
	dQ = [0.0]*len(Q-1)

	for i in xrange(1,len(Q)):
		if(i==0):
			dQ[i] = Q[i] + ( (I_vol[i]/V0)*((Q[i]-0.0000)/2) )
		else:
			dQ[i] = Q[i] + ( (I_vol[i]/V0)*((Q[i]+Q[i-1])/2) ) - Q[i-1]

	return dQ