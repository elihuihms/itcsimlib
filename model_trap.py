from thermo		import *
from itc_model	import ITCModel
from itc_calc	import ITCCalcDLL

class NN(ITCModel):
	"""
	A twelve-parameter model describing all possible nearest-neighbor interactions
	"""

	def __init__(self):
		ITCModel.__init__(self)

		self.add_param( 'dG0',	'dG',	description='Intrinsic free energy change upon binding' )
		self.add_param( 'dGa',	'dG',	description='Free energy of coupling to two unoccupied sites' )
		self.add_param( 'dGb',	'dG',	description='Free energy of coupling to one unoccupied site' )
		self.add_param( 'dGc',	'dG',	description='Free energy of coupling to two occupied sites' )
		self.add_param( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding' )
		self.add_param( 'dHa',	'dH',	description='Enthalpy of coupling to two unoccupied sites' )
		self.add_param( 'dHb',	'dH',	description='Enthalpy of coupling to one unoccupied site' )
		self.add_param( 'dHc',	'dH',	description='Enthalpy of coupling to two occupied sites' )
		self.add_param( 'dCp0',	'dCp',	description='Intrinsic change in heat capacity change upon binding' )
		self.add_param( 'dCpa',	'dCp',	description='Change in heat capacity of coupling to two unoccupied sites' )
		self.add_param( 'dCpb',	'dCp',	description='Change in heat capacity of coupling to one unoccupied site' )
		self.add_param( 'dCpc',	'dCp',	description='Change in heat capacity of coupling to two occupied sites' )

	def get_format(self, T, T_ref):
		return [
			dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T_ref ),
			dG_vant_Hoff( self.params['dGa'], self.params['dHa'], self.params['dCpa'], T, T_ref ),
			dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T_ref ),
			dG_vant_Hoff( self.params['dGc'], self.params['dHc'], self.params['dCpc'], T, T_ref ),
			dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T_ref ),
			dH_vant_Hoff( self.params['dHa'], self.params['dCpa'], T, T_ref ),
			dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T_ref ),
			dH_vant_Hoff( self.params['dHc'], self.params['dCpc'], T, T_ref )
		]

	def get_worker(self,in_queue,out_queue):
		return ITCCalcDLL('./itcsimlib/model_trap_nn.so',[11,1],in_queue,out_queue)

class SK(ITCModel):
	"""
	A nine-parameter model describing the Saroff-Kiefer model
	"""

	def __init__(self):
		ITCModel.__init__(self)

		self.add_param( 'dG0',	'dG',	description='Intrinsic free energy change upon binding' )
		self.add_param( 'dGa',	'dG',	description='Free energy of coupling to an unoccupied site' )
		self.add_param( 'dGb',	'dG',	description='Free energy of coupling to an occupied site' )
		self.add_param( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding' )
		self.add_param( 'dHa',	'dH',	description='Enthalpy of coupling to an unoccupied site' )
		self.add_param( 'dHb',	'dH',	description='Enthalpy of coupling to an occupied site' )
		self.add_param( 'dCp0',	'dCp',	description='Intrinsic change in heat capacity change upon binding' )
		self.add_param( 'dCpa',	'dCp',	description='Change in heat capacity of coupling to an unoccupied site' )
		self.add_param( 'dCpb',	'dCp',	description='Change in heat capacity of coupling to an occupied site' )

	def get_format(self, T, T_ref):
		return [
			dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T_ref ),
			dG_vant_Hoff( self.params['dGa'], self.params['dHa'], self.params['dCpa'], T, T_ref ),
			dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T_ref ),
			dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T_ref ),
			dH_vant_Hoff( self.params['dHa'], self.params['dCpa'], T, T_ref ),
			dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T_ref ),
		]

	def get_worker(self,in_queue,out_queue):
		return ITCCalcDLL('./itcsimlib/model_trap_sk.so',[11,1],in_queue,out_queue)

class IK(ITCModel):
	"""
	A nine-parameter model describing the Kleckner model
	"""

	def __init__(self):
		ITCModel.__init__(self)

		self.add_param( 'dGX',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_param( 'dGY',	'dG',	description='Free energy change upon binding to a site flanked by one occupied' )
		self.add_param( 'dGZ',	'dG',	description='Free energy change upon binding to a site flanked by two occupied' )
		self.add_param( 'dHX',	'dH',	description='Enthalpy change upon binding to a site flanked by two unoccupied' )
		self.add_param( 'dHY',	'dH',	description='Enthalpy change upon binding to a site flanked by one occupied' )
		self.add_param( 'dHZ',	'dH',	description='Enthalpy change upon binding to a site flanked by two occupied' )
		self.add_param( 'dCpX',	'dCp',	description='Change in heat capacity of coupling to a site flanked by two unoccupied' )
		self.add_param( 'dCpY',	'dCp',	description='Change in heat capacity of coupling to a site flanked by one occupied' )
		self.add_param( 'dCpZ',	'dCp',	description='Change in heat capacity of coupling to a site flanked by two occupied' )

	def get_format(self, T, T_ref):
		return [
			dG_vant_Hoff( self.params['dGX'], self.params['dHX'], self.params['dCpX'], T, T_ref ),
			dG_vant_Hoff( self.params['dGY'], self.params['dHY'], self.params['dCpY'], T, T_ref ),
			dG_vant_Hoff( self.params['dGZ'], self.params['dHZ'], self.params['dCpZ'], T, T_ref ),
			dH_vant_Hoff( self.params['dHX'], self.params['dCpX'], T, T_ref ),
			dH_vant_Hoff( self.params['dHY'], self.params['dCpY'], T, T_ref ),
			dH_vant_Hoff( self.params['dHZ'], self.params['dCpZ'], T, T_ref ),
		]

	def get_worker(self,in_queue,out_queue):
		return ITCCalcDLL('./itcsimlib/model_trap_ik.so',[11,1],in_queue,out_queue)
	
class IKi(ITCModel):
	"""
	A nine-parameter model describing the Kleckner model with an intrinsic binding coeff
	"""

	def __init__(self):
		ITCModel.__init__(self)

		self.add_param( 'dGi',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_param( 'dGoe',	'dG',	description='Additional free energy change upon binding to a site flanked by one occupied' )
		self.add_param( 'dGoo',	'dG',	description='Additional free energy change upon binding to a site flanked by two occupied' )
		self.add_param( 'dHi',	'dH',	description='Enthalpy change upon binding to a site flanked by two unoccupied' )
		self.add_param( 'dHoe',	'dH',	description='Additional enthalpy change upon binding to a site flanked by one occupied' )
		self.add_param( 'dHoo',	'dH',	description='Additional enthalpy change upon binding to a site flanked by two occupied' )
		self.add_param( 'dCpi',	'dCp',	description='Change in heat capacity of coupling to a site flanked by two unoccupied' )
		self.add_param( 'dCpoe',	'dCp',	description='Additional change in heat capacity of coupling to a site flanked by one occupied' )
		self.add_param( 'dCpoo',	'dCp',	description='Additional change in heat capacity of coupling to a site flanked by two occupied' )

	def get_format(self, T, T_ref):
		return [
			dG_vant_Hoff( self.params['dGi'], self.params['dHi'], self.params['dCpi'], T, T_ref ),
			dG_vant_Hoff( self.params['dGi']+self.params['dGoe'], self.params['dHi']+self.params['dHoe'], self.params['dCpi']+self.params['dCpoe'], T, T_ref ),
			dG_vant_Hoff( self.params['dGi']+self.params['dGoo'], self.params['dHi']+self.params['dHoo'], self.params['dCpi']+self.params['dCpoo'], T, T_ref ),
			dH_vant_Hoff( self.params['dHi'], self.params['dCpi'], T, T_ref ),
			dH_vant_Hoff( self.params['dHi']+self.params['dHoe'], self.params['dCpi']+self.params['dCpoe'], T, T_ref ),
			dH_vant_Hoff( self.params['dHi']+self.params['dHoo'], self.params['dCpi']+self.params['dCpoo'], T, T_ref ),
		]

	def get_worker(self,in_queue,out_queue):
		return ITCCalcDLL('./itcsimlib/model_trap_ik.so',[11,1],in_queue,out_queue)