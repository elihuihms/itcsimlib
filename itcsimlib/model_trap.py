"""Binding models for the undecameric TRAP + tryptophan system, demonstrating the use of models written in C via DLLs.

Using compiled DLLs typically provides an order-of-magnitude increase or better in execution speed over native Python models.

"""

import os
import scipy
import ctypes

from .itc_model	import ITCModel
from .thermo	import *


class TRAP_DLL_Model(ITCModel):
	"""A model that uses a shared object library to calculate evolved heats for the TRAP + Tryptophan system.

	Attributes
	----------
	libpath : string
		The path to the dll, relative to the directory this module is present in. 
	"""

	libpath = None

	def __init__(self):
		ITCModel.__init__(self)
		
		self.nsites,self.circular = 11,1

		self._lib = None
		self._path = os.path.join( os.path.dirname(__file__), self.libpath )

		if not os.path.exists(self._path):
			raise ImportError("Could not import shared library path at \"%s\"."%self._path)

		# Dry run to attempt to load the library.
		self.start()
		self.stop()

		self.add_component('TRAP',description='An %i-site circular lattice of tryptophan binding sites'%(self.nsites))
		self.add_component('Trp',description='A molecule of tryptophan')

	def start(self):
		"""Loads the specified shared library.

		Returns
		-------
		errno : integer
			The return code from the DLL, 0 if no error.
		"""
		self._lib = ctypes.cdll.LoadLibrary(self._path)
		return self._lib.setup(ctypes.c_int(self.nsites),ctypes.c_int(self.circular))
		
	def stop(self):
		"""Closes the shared library."""
		return self._lib.close()

	def calc(self,T,concentrations,params):
		"""Converts the arguments to their appropriate ctypes, passes to the DLL, and returns the calculated enthalpies.

		Arguments
		---------
		T : float
			The experimental temperature (in Kelvin)
		concentrations : list of dicts
			Concentrations of components from the experiment at each injection point.
		params : list of floats
			Model parameter values

		Returns
		-------
		list of doubles
			The integrated enthalpies at each injection point.

		Raises
		------
		Exception
			If DLL returns a non-zero error code.

		Notes
		-----
			Named components are "TRAP" (lattice) and "Trp" (ligand).
		"""
		n = len(concentrations)
		Q = scipy.zeros(n,scipy.dtype('d'))
		
		# patch for compatibility with general model nomenclature
		if 'TRAP' not in concentrations[0]:
			concentrations = [{'TRAP':concentrations[i]['Macromolecule'],'Trp':concentrations[i]['Ligand']} for i in range(n)]
				
		status = self._lib.calc(
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
	A nine-parameter model describing the additive 1997 Saroff-Kiefer model, as published.

	Notes
	-----
		The parameters dG, dGa, and dBa (and their enthalpic and heat capacity counterparts) are correlated and should not be fit simultaneously, as described in Ihms. et al., Biophysical Journal (2017), 
	"""

	libpath = 'model_trap_sk.so'

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
		
	def Q(self,T0,T,concentrations):
		"""Returns the total binding heat at each injection predicted by the model and its current parameter values. See parent model for information."""
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
	A nine-parameter model describing the generalized Kleckner (non-additive) model of binding.
	"""

	libpath = 'model_trap_ik.so'

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
		
	def Q(self,T0,T,concentrations):
		"""Returns the total binding heat at each injection predicted by the model and its current parameter values. See parent model for information."""
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
	A nine-parameter model describing the Kleckner (non-additive) model, expressed using an intrinsic binding coeffcient.
	"""

	libpath = 'model_trap_ik.so'

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

	def Q(self,T0,T,concentrations):
		"""Returns the total binding heat at each injection predicted by the model and its current parameter values. See parent model for information."""
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
	A six-parameter model describing the additive 1997 Saroff-Kiefer model, using an intrinsic binding coefficient to remove correlated terms.
	"""

	libpath = 'model_trap_ik.so'

	def __init__(self):
		TRAP_DLL_Model.__init__(self)

		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding' )
		self.add_parameter( 'dGb',	'dG',	description='Free energy of coupling to an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding' )
		self.add_parameter( 'dHb',	'dH',	description='Enthalpy of coupling to an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic change in heat capacity change upon binding' )
		self.add_parameter( 'dCpb',	'dCp',	description='Change in heat capacity of coupling to an occupied site' )

	def Q(self,T0,T,concentrations):
		"""Returns the total binding heat at each injection predicted by the model and its current parameter values. See parent model for information."""
		p = (
			dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 ),
			dG_vant_Hoff( self.params['dG0']+(1*self.params['dGb']), self.params['dH0']+(1*self.params['dHb']), self.params['dCp0']+(1*self.params['dCpb']), T, T0 ),
			dG_vant_Hoff( self.params['dG0']+(2*self.params['dGb']), self.params['dH0']+(2*self.params['dHb']), self.params['dCp0']+(2*self.params['dCpb']), T, T0 ),
			dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 ),
			dH_vant_Hoff( self.params['dH0']+(1*self.params['dHb']), self.params['dCp0']+(1*self.params['dCpb']), T, T0 ),
			dH_vant_Hoff( self.params['dH0']+(2*self.params['dHb']), self.params['dCp0']+(2*self.params['dCpb']), T, T0 ),
			)
		return self.calc(T,concentrations,p)
