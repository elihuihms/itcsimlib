from math		import exp,pow,sqrt
import scipy.optimize as optimize

from itc_model	import ITCModel
from thermo		import *

class OneMode(ITCModel):
	"""
	A four-parameter phenomological model describing binding to a single site type
	"""

	def __init__(self):
		ITCModel.__init__(self)

		self.add_component('Macromolecule')
		self.add_component('Ligand')

		self.add_parameter( 'n',	'n',	description='n_sites', bounds=[0,None], default=1.0 )
		self.add_parameter( 'dG',	'dG',	description='Free energy change upon binding' )
		self.add_parameter( 'dH',	'dH',	description='Enthalpy change upon binding' )
		self.add_parameter( 'dCp',	'dCp',	description='Heat capacity change' )

	def Q(self,T0,T,concentrations):	
		n1,Ka,dH = (
			self.params['n'],
			1.0/Kd_from_dG( dG_vant_Hoff( self.params['dG'], self.params['dH'], self.params['dCp'], T, T0 ), T),
			dH_vant_Hoff( self.params['dH'], self.params['dCp'], T, T0 )
		)

		Q = [0.0]*len(concentrations)
		for i,c in enumerate(concentrations):
			Q[i] = ((n1*dH)/2.0)*(1.0 +(c['Ligand']/(n1*c['Macromolecule'])) +(1.0/(n1*Ka*c['Macromolecule'])) -sqrt( pow(1 +(c['Ligand']/(n1*c['Macromolecule'])) +(1.0/(n1*Ka*c['Macromolecule'])), 2.0) -((4.0*c['Ligand'])/(n1*c['Macromolecule']))))
		return Q

class NModes(ITCModel):
	"""
	A 4n-parameter phenomological model describing binding to n independent types of sites
	"""

	def __init__(self,nmodes=2):
		ITCModel.__init__(self)
		self.nmodes = nmodes

		self.add_component('Macromolecule')
		self.add_component('Ligand')

		for i in xrange(self.nmodes):
			self.add_parameter( "n%i"%(i),		'n',	description='n_sites', bounds=[0,None], default=1.0 )
			self.add_parameter( "dG%i"%(i),		'dG',	description='Free energy change upon binding' )
			self.add_parameter( "dH%i"%(i),		'dH',	description='Enthalpy change upon binding' )
			self.add_parameter( "dCp%i"%(i),	'dCp',	description='Heat capacity change' )

	def Q(self,T0,T,concentrations):
		n,p = len(concentrations),[None]*(self.nmodes*3)
		for i in xrange(self.nmodes):
			n,dG,dH,dCp	= 'n'+str(i),'dG'+str(i),'dH'+str(i),'dCp'+str(i)
			p[i*3 +0]	= self.params[n]
			p[i*3 +1]	= 1.0/Kd_from_dG( dG_vant_Hoff( self.params[dG], self.params[dH], self.params[dCp], T, T0 ), T)
			p[i*3 +2]	= dH_vant_Hoff( self.params[dH], self.params[dCp], T, T0 )

		def _get_free(Lfree,Ltot,Ptot):
			Lbound = 0.0
			for i in xrange(self.nmodes):
				n,Ka,dH = p[i*3:i*3+3]
				Lbound += n * Ptot * (Ka*Lfree)/(Ka*Lfree +1)
			return Ltot -Lbound -Lfree

		Q = [0.0]*n
		for j,c in enumerate(concentrations):
			Lfree = optimize.brentq( _get_free, 0.0, c['Ligand'], args=(c['Ligand'],c['Macromolecule'],p) )
			for i in xrange(len(self.params)/3):
				n,Ka,dH = self.params[i*3:i*3+3]
				Q[j] +=( n * dH * (Ka*Lfree)/(Ka*Lfree +1) )
		return Q
