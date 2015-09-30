from math		import exp,pow,sqrt

from thermo		import *
from itc_model	import ITCModel
from itc_calc	import ITCCalc

from ctypes		import c_double

class OneMode(ITCModel):
	"""
	A four-parameter phenomological model describing binding to a single site type
	"""

	def __init__(self):
		ITCModel.__init__(self)

		self.add_param( 'n',	'n',	description='n_sites', bounds=[0,None], default=1.0 )
		self.add_param( 'dG',	'dG',	description='Free energy change upon binding' )
		self.add_param( 'dH',	'dH',	description='Enthalpy change upon binding' )
		self.add_param( 'dCp',	'dCp',	description='Heat capacity change' )

	def get_params(self, T, T_ref):
		return [
			self.params['n'],
			1.0/Kd_from_dG( dG_vant_Hoff( self.params['dG'], self.params['dH'], self.params['dCp'], T, T_ref ), T),
			dH_vant_Hoff( self.params['dH'], self.params['dCp'], T, T_ref )
		]

	def get_worker(self,in_queue,out_queue):
		def _calc(x,T,P,L,Q,params):
			n1,Ka,dH = params

			""" binding polynomial from Microcal's "ITC Data Analysis In Origin, Rev. F-4, p.g. 106, eq. 9 """
			for i in xrange(x):
				Q[i] = ((n1*dH)/2.0)*(1.0 +(L[i]/(n1*P[i])) +(1.0/(n1*Ka*P[i])) -sqrt( pow(1 +(L[i]/(n1*P[i])) +(1.0/(n1*Ka*P[i])), 2.0) -((4.0*L[i])/(n1*P[i]))))

		return ITCCalc(_calc,in_queue,out_queue)

class NModes(ITCModel):
	"""
	A 4n-parameter phenomological model describing binding to n independent types of sites
	"""

	def __init__(self,nmodes=2):
		ITCModel.__init__(self)
		self.nmodes = nmodes

		for i in xrange(self.nmodes):
			self.add_param( "n%i"%(i),		'n',	description='n_sites', bounds=[0,None], default=1.0 )
			self.add_param( "dG%i"%(i),		'dG',	description='Free energy change upon binding' )
			self.add_param( "dH%i"%(i),		'dH',	description='Enthalpy change upon binding' )
			self.add_param( "dCp%i"%(i),	'dCp',	description='Heat capacity change' )

	def get_format(self, T, T_ref):
		ret = [None]*(self.nmodes*3)
		for i in xrange(self.nmodes):
			n,dG,dH,dCp	= 'n'+str(i),'dG'+str(i),'dH'+str(i),'dCp'+str(i)
			ret[i*3 +0]	= self.params[n]
			ret[i*3 +1]	= 1.0/Kd_from_dG( dG_vant_Hoff( self.params[dG], self.params[dH], self.params[dCp], T, T_ref ), T)
			ret[i*3 +2]	= dH_vant_Hoff( self.params[dH], self.params[dCp], T, T_ref )
		return ret

	def get_worker(self,in_queue,out_queue):
		def _calc(x,T,P,L,Q,params):
			import scipy.optimize as optimize

			def _get_free(Lfree,Ltot,Ptot,params):
				Lbound = 0
				for i in xrange(len(params)/3):
					n,Ka,dH = params[i*3:i*3+3]
					Lbound += n * Ptot * (Ka*Lfree)/(Ka*Lfree +1)

				return Ltot -Lbound -Lfree

			for j in xrange(x):
				Lfree = optimize.brentq( _get_free, 0, L[j], args=(L[j],P[j],params) )
				for i in xrange(len(params)/3):
					n,Ka,dH = params[i*3:i*3+3]
					Q[j] +=( n * dH * (Ka*Lfree)/(Ka*Lfree +1) )

		return ITCCalc(_calc,in_queue,out_queue)

#class OneModeC(ITCModel):
#	"""
#	A four-parameter phenomological model describing binding to a single site type
#	"""
#
#	def __init__(self):
#		ITCModel.__init__(self)
#
#		self.add_param( 'n',	'n',	description='n_sites', bounds=[0,None], default=1.0 )
#		self.add_param( 'dG',	'dG',	description='Free energy change upon binding' )
#		self.add_param( 'dH',	'dH',	description='Enthalpy change upon binding' )
#		self.add_param( 'dCp',	'dCp',	description='Heat capacity change' )
#
#	def get_params(self, T, T_ref):
#		return [
#			c_double(self.get_param('n')),
#			c_double(dG_vant_Hoff( self.get_param('dG'), self.get_param('dH'), self.get_param('dCp'), T, T_ref )),
#			c_double(dH_vant_Hoff( self.get_param('dH'), self.get_param('dCp'), T, T_ref ))
#		]
#
#	def get_worker(self,in_queue,out_queue):
#		return ITCCalcDLL('./itcsimlib/itc_sim_n_identical.so',[],in_queue,out_queue)



