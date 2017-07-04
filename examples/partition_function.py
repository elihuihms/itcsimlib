import sympy
sympy.init_printing()

from itcsimlib import *
from itcsimlib.model_ising import *

print sympy.latex(NonAdditive(nsites=11,circular=True).get_partition_function())
print

class Hemoglobin4P(Ising):
	def __init__(self,nsites=4,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG1',	'dG',	description='Free energy of binding the first O2 molecule to hemoglobin' )
		self.add_parameter( 'dG2',	'dG',	description='Free energy of binding the second O2 molecule to hemoglobin' )
		self.add_parameter( 'dG3',	'dG',	description='Free energy of binding the third O2 molecule to hemoglobin' )
		self.add_parameter( 'dG4',	'dG',	description='Free energy of binding the fourth O2 molecule to hemoglobin' )
		self.add_parameter( 'dH1',	'dH',	description='Enthalpy of binding the first O2 molecule to hemoglobin' )
		self.add_parameter( 'dH2',	'dH',	description='Enthalpy of binding the second O2 molecule to hemoglobin' )
		self.add_parameter( 'dH3',	'dH',	description='Enthalpy of binding the third O2 molecule to hemoglobin' )
		self.add_parameter( 'dH4',	'dH',	description='Enthalpy of binding the fourth O2 molecule to hemoglobin' )
		self.add_parameter( 'dCp1',	'dCp',	description='Heat capacity change of binding the first O2 molecule to hemoglobin' )
		self.add_parameter( 'dCp2',	'dCp',	description='Heat capacity change of binding the second O2 molecule to hemoglobin' )
		self.add_parameter( 'dCp3',	'dCp',	description='Heat capacity change of binding the third O2 molecule to hemoglobin' )
		self.add_parameter( 'dCp4',	'dCp',	description='Heat capacity change of binding the fourth O2 molecule to hemoglobin' )

	def set_energies(self,T0,T):
		dG1 = dG_vant_Hoff( self.params['dG1'], self.params['dH1'], self.params['dCp1'], T, T0 )
		dG2 = dG_vant_Hoff( self.params['dG2'], self.params['dH2'], self.params['dCp2'], T, T0 )
		dG3 = dG_vant_Hoff( self.params['dG3'], self.params['dH3'], self.params['dCp3'], T, T0 )
		dG4 = dG_vant_Hoff( self.params['dG4'], self.params['dH4'], self.params['dCp4'], T, T0 )
		dH1 = dH_vant_Hoff( self.params['dH1'], self.params['dCp1'], T, T0 )
		dH2 = dH_vant_Hoff( self.params['dH2'], self.params['dCp2'], T, T0 )
		dH3 = dH_vant_Hoff( self.params['dH3'], self.params['dCp3'], T, T0 )
		dH4 = dH_vant_Hoff( self.params['dH4'], self.params['dCp4'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			bound = 0 # for this config, how many bound ligands?
			for j in xrange(self.nsites):
				if self.get_site_occupancy(i,j):
					bound +=1
			
			if bound == 1:
				self.gibbs[i] = dG1
				self.enthalpies[i] = dH1
				self.config_expressions[i] += self.parameter_symbols['dG1']
			elif bound == 2:
				self.gibbs[i] = dG2
				self.enthalpies[i] = dH2
				self.config_expressions[i] += self.parameter_symbols['dG2']
			elif bound == 3:
				self.gibbs[i] = dG3
				self.enthalpies[i] = dH3
				self.config_expressions[i] += self.parameter_symbols['dG3']
			elif bound == 4:
				self.gibbs[i] = dG4
				self.enthalpies[i] = dH4
				self.config_expressions[i] += self.parameter_symbols['dG4']

print sympy.latex(Hemoglobin4P().get_partition_function())
print

class Hemoglobin3P(Ising):
	def __init__(self,nsites=4,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy of binding' )
		self.add_parameter( 'dGa',	'dG',	description='Free energy of coupling to an empty site' )
		self.add_parameter( 'dGb',	'dG',	description='Free energy of coupling to an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy of binding' )
		self.add_parameter( 'dHa',	'dH',	description='Enthalpy of coupling to an empty site' )
		self.add_parameter( 'dHb',	'dH',	description='Enthalpy of coupling to an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic heat capacity change of binding' )
		self.add_parameter( 'dCpa',	'dCp',	description='Heat capacity change of coupling to an empty site' )
		self.add_parameter( 'dCpb',	'dCp',	description='Heat capacity change of coupling to an occupied site' )

	def set_energies(self,T0,T):
		dG0 = dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 )
		dGa = dG_vant_Hoff( self.params['dGa'], self.params['dHa'], self.params['dCpa'], T, T0 )
		dGb = dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 )
		dH0 = dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 )
		dHa = dH_vant_Hoff( self.params['dHa'], self.params['dCpa'], T, T0 )
		dHb = dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			for j in xrange(self.nsites):

				if self.get_site_occupancy(i,j): # is the site occupied?
					self.gibbs[i] += dG0
					self.enthalpies[i] += dH0
					self.config_expressions[i] += self.parameter_symbols['dG0']

					if self.get_site_occupancy(i,j+1): # is the next site occupied?
						self.gibbs[i] += dGb
						self.enthalpies[i] += dHb
						self.config_expressions[i] += self.parameter_symbols['dGb']
					else:
						self.gibbs[i] += dGa
						self.enthalpies[i] += dHa
						self.config_expressions[i] += self.parameter_symbols['dGa']

					if self.get_site_occupancy(i,j-1): # is the previous site occupied?
						self.gibbs[i] += dGb
						self.enthalpies[i] += dHb
						self.config_expressions[i] += self.parameter_symbols['dGb']
					else:
						self.gibbs[i] += dGa
						self.enthalpies[i] += dHa
						self.config_expressions[i] += self.parameter_symbols['dGa']

print sympy.latex(Hemoglobin3P().get_partition_function())
print

class Hemoglobin2P(Ising):
	def __init__(self,nsites=4,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy of binding' )
		self.add_parameter( 'dGb',	'dG',	description='Free energy of coupling to an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy of binding' )
		self.add_parameter( 'dHb',	'dH',	description='Enthalpy of coupling to an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic heat capacity change of binding' )
		self.add_parameter( 'dCpb',	'dCp',	description='Heat capacity change of coupling to an occupied site' )

	def set_energies(self,T0,T):
		self.precision = 1E-20
		dG0 = dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 )
		dGb = dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 )
		dH0 = dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 )
		dHb = dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			for j in xrange(self.nsites):

				if self.get_site_occupancy(i,j): # is the site occupied?
					self.gibbs[i] += dG0
					self.enthalpies[i] += dH0
					self.config_expressions[i] += self.parameter_symbols['dG0']

					if self.get_site_occupancy(i,j+1): # is the next site occupied?
						self.gibbs[i] += dGb
						self.enthalpies[i] += dHb
						self.config_expressions[i] += self.parameter_symbols['dGb']

					if self.get_site_occupancy(i,j-1): # is the previous site occupied?
						self.gibbs[i] += dGb
						self.enthalpies[i] += dHb
						self.config_expressions[i] += self.parameter_symbols['dGb']

print sympy.latex(Hemoglobin2P().get_partition_function())
print

sim = ITCSim(T0=273.15+25,verbose=True,threads=1,units='kcal')
sim.add_experiment_synthetic(
	T=298.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Lattice":1E-6},
	Syringe={"Ligand":40E-6},
	title='Hemoglobin2P_example')

sim.set_model(Hemoglobin2P())
sim.set_model_params(dG0=-8,dGb=-0.5,dH0=-10,dHb=-2)
sim.run()
sim.make_plots(hardcopy=True)
sim.done()