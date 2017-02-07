import os
import scipy
import ctypes

from itc_model	import ITCModel
from thermo		import *

class TRAP_DLL_Model(ITCModel):
	
	def __init__(self):
		ITCModel.__init__(self)
		
		self.nsites,self.circular = 11,1
		
		self.add_component('TRAP',description='An %i-site circular lattice of tryptophan binding sites'%(self.nsites))
		self.add_component('Trp',description='A molecule of tryptophan')
	
	def loadDLL(self,lib):
		self.lib = ctypes.cdll.LoadLibrary(os.path.join( os.path.dirname(__file__), lib ))
		return self.lib.setup(ctypes.c_int(self.nsites),ctypes.c_int(self.circular))
		
	def stop(self):
		return self.lib.close()

	def calc(self,T,concentrations,params):
		n = len(concentrations)
		Q = scipy.zeros(n,scipy.dtype('d'))
		
		# patch for compatibility with general model nomenclature
		if 'TRAP' not in concentrations[0]:
			concentrations = [{'TRAP':concentrations[i]['Macromolecule'],'Trp':concentrations[i]['Ligand']} for i in xrange(n)]
				
		status = self.lib.calc(
			ctypes.c_int( n ),
			ctypes.c_double( T ),
			scipy.array([c['TRAP'] for c in concentrations],scipy.dtype('d')).ctypes,
			scipy.array([c['Trp']  for c in concentrations],scipy.dtype('d')).ctypes,
			Q.ctypes,
			scipy.array(params,scipy.dtype('d')).ctypes
		)
		
		if status != 0:
			raise Exception("DLL returned a non-zero error code: %i"%(status))

		return Q

class SK(TRAP_DLL_Model):
	"""
	A nine-parameter model describing the Saroff-Kiefer model
	"""

	def __init__(self):
		TRAP_DLL_Model.__init__(self)

		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding' )
		self.add_parameter( 'dGa',	'dG',	description='Free energy of coupling to an unoccupied site' )
		self.add_parameter( 'dGb',	'dG',	description='Free energy of coupling to an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding' )
		self.add_parameter( 'dHa',	'dH',	description='Enthalpy of coupling to an unoccupied site' )
		self.add_parameter( 'dHb',	'dH',	description='Enthalpy of coupling to an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic change in heat capacity change upon binding' )
		self.add_parameter( 'dCpa',	'dCp',	description='Change in heat capacity of coupling to an unoccupied site' )
		self.add_parameter( 'dCpb',	'dCp',	description='Change in heat capacity of coupling to an occupied site' )

	def start(self):
		return self.loadDLL('model_trap_sk.so')
		
	def Q(self,T0,T,concentrations):
		p = (
			dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 ),
			dG_vant_Hoff( self.params['dGa'], self.params['dHa'], self.params['dCpa'], T, T0 ),
			dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 ),
			dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 ),
			dH_vant_Hoff( self.params['dHa'], self.params['dCpa'], T, T0 ),
			dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
			)
		return self.calc(T,concentrations,p)

class IK(TRAP_DLL_Model):
	"""
	A nine-parameter model describing the Kleckner model
	"""

	def __init__(self):
		TRAP_DLL_Model.__init__(self)

		self.add_parameter( 'dGX',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dGY',	'dG',	description='Free energy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dGZ',	'dG',	description='Free energy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dHX',	'dH',	description='Enthalpy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dHY',	'dH',	description='Enthalpy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dHZ',	'dH',	description='Enthalpy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dCpX',	'dCp',	description='Change in heat capacity of coupling to a site flanked by two unoccupied' )
		self.add_parameter( 'dCpY',	'dCp',	description='Change in heat capacity of coupling to a site flanked by one occupied' )
		self.add_parameter( 'dCpZ',	'dCp',	description='Change in heat capacity of coupling to a site flanked by two occupied' )

	def start(self):
		return self.loadDLL('model_trap_ik.so')
		
	def Q(self,T0,T,concentrations):
		p = (
			dG_vant_Hoff( self.params['dGX'], self.params['dHX'], self.params['dCpX'], T, T0 ),
			dG_vant_Hoff( self.params['dGY'], self.params['dHY'], self.params['dCpY'], T, T0 ),
			dG_vant_Hoff( self.params['dGZ'], self.params['dHZ'], self.params['dCpZ'], T, T0 ),
			dH_vant_Hoff( self.params['dHX'], self.params['dCpX'], T, T0 ),
			dH_vant_Hoff( self.params['dHY'], self.params['dCpY'], T, T0 ),
			dH_vant_Hoff( self.params['dHZ'], self.params['dCpZ'], T, T0 )
			)
		return self.calc(T,concentrations,p)
	
class IKi(TRAP_DLL_Model):
	"""
	A nine-parameter model describing the Kleckner model with an intrinsic binding coeff
	"""

	def __init__(self):
		TRAP_DLL_Model.__init__(self)

		self.add_parameter( 'dG0',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dGoe',	'dG',	description='Additional free energy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dGoo',	'dG',	description='Additional free energy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dH0',	'dH',	description='Enthalpy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dHoe',	'dH',	description='Additional enthalpy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dHoo',	'dH',	description='Additional enthalpy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dCp0',	'dCp',	description='Change in heat capacity of coupling to a site flanked by two unoccupied' )
		self.add_parameter( 'dCpoe','dCp',	description='Additional change in heat capacity of coupling to a site flanked by one occupied' )
		self.add_parameter( 'dCpoo','dCp',	description='Additional change in heat capacity of coupling to a site flanked by two occupied' )

	def start(self):
		return self.loadDLL('model_trap_ik.so')

	def Q(self,T0,T,concentrations):
		p = (
			dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 ),
			dG_vant_Hoff( self.params['dG0']+self.params['dGoe'], self.params['dH0']+self.params['dHoe'], self.params['dCp0']+self.params['dCpoe'], T, T0 ),
			dG_vant_Hoff( self.params['dG0']+self.params['dGoo'], self.params['dH0']+self.params['dHoo'], self.params['dCp0']+self.params['dCpoo'], T, T0 ),
			dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 ),
			dH_vant_Hoff( self.params['dH0']+self.params['dHoe'], self.params['dCp0']+self.params['dCpoe'], T, T0 ),
			dH_vant_Hoff( self.params['dH0']+self.params['dHoo'], self.params['dCp0']+self.params['dCpoo'], T, T0 )
			)
		return self.calc(T,concentrations,p)
		
class SKa(TRAP_DLL_Model):
	"""
	A Six-parameter model describing an additive model of cooperativity 
	"""

	def __init__(self):
		TRAP_DLL_Model.__init__(self)

		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding' )
		self.add_parameter( 'dGb',	'dG',	description='Free energy of coupling to an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding' )
		self.add_parameter( 'dHb',	'dH',	description='Enthalpy of coupling to an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic change in heat capacity change upon binding' )
		self.add_parameter( 'dCpb',	'dCp',	description='Change in heat capacity of coupling to an occupied site' )

	def start(self):
		return self.loadDLL('model_trap_ik.so')

	def Q(self,T0,T,concentrations):
		p = (
			dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 ),
			dG_vant_Hoff( self.params['dG0']+(1*self.params['dGb']), self.params['dH0']+(1*self.params['dHb']), self.params['dCp0']+(1*self.params['dCpb']), T, T0 ),
			dG_vant_Hoff( self.params['dG0']+(2*self.params['dGb']), self.params['dH0']+(2*self.params['dHb']), self.params['dCp0']+(2*self.params['dCpb']), T, T0 ),
			dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 ),
			dH_vant_Hoff( self.params['dH0']+(1*self.params['dHb']), self.params['dCp0']+(1*self.params['dCpb']), T, T0 ),
			dH_vant_Hoff( self.params['dH0']+(2*self.params['dHb']), self.params['dCp0']+(2*self.params['dCpb']), T, T0 ),
			)
		return self.calc(T,concentrations,p)
