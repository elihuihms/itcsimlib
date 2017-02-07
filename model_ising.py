import math
import scipy
import scipy.optimize

from itc_model	import ITCModel
from thermo		import _R
from thermo		import *

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

class Ising(ITCModel):

	def __init__(self,nsites=3,circular=1):
		ITCModel.__init__(self)
		
		self.nsites,self.circular = nsites,circular
		
		self.nconfigs	= 2**self.nsites
		self.configs	= [ [int(s) for s in ("{0:0%ib}"%(self.nsites)).format(i)] for i in xrange(self.nconfigs) ]
		self.config_params	= [ [] for i in xrange(self.nconfigs) ] # list of free energy parameters used to construct the model's partition function
		self.bound		= [ c.count(1) for c in self.configs ] # number of bound sites
		self.weights	= [0.0]*self.nconfigs # probability of each config
		self.gibbs		= [0.0]*self.nconfigs # free energy of each config
		self.enthalpies	= [0.0]*self.nconfigs # enthalpic energy of each config
		self.precision	= 1E-9 # the precision in ligand concentration required for convergence during set_probabilities()

		if self.circular:		
			self.add_component('Lattice',description='A circular lattice with %i binding sites'%(self.nsites))
		else:
			self.add_component('Lattice',description='A linear lattice with %i binding sites'%(self.nsites))
		self.add_component('Ligand',description='A lattice-binding ligand')
		
		return

	def get_site_occupancy(self,config,site):
		if site < 0:
			if not self.circular:
				return None
			return self.configs[config][site+self.nsites] == 1
		if site >= self.nsites:
			if not self.circular:
				return None
			return self.configs[config][site%self.nsites] == 1
		return self.configs[config][site] == 1
		
	def get_partition_function(self):
		self.set_energies(273.15,273.15)
		
		ret = "\\begin{array}{lcl}\n"
		ret+= "\\Xi\n"
		ret+= "&=&1\\\\\n"
		for n in xrange(1,self.nsites+1):
			if n == 1:
				ret+= "&+& L*"
			elif n > 1:
				ret+= "&+& L^{%i}*"%(n)
				
			config_terms = [ sorted(self.config_params[i]) for i in xrange(self.nconfigs) if self.bound[i] == n]
			config_gibbs = [ self.gibbs[i] for i in xrange(self.nconfigs) if self.bound[i] == n ]
			
			# figure out what parameters are common in all terms for these config expressions
			# add them to the partition function, and then remove them from the config lists
			for param in set([parameter for terms in config_terms for parameter in terms]):
				num = min([term.count(param) for term in config_terms])
				for j in xrange(num):
					for i in xrange(len(config_terms)):				
						config_terms[i].remove(param)
				
				if num == 1:
					ret+= "%s*"%(param)
				elif num > 1:
					ret+= "%s^{%i}*"%(param,num)

			ret = ret[:-1] # remove last *
			ret+= "["
			
			# condense degenerate configuration expressions
			param_sets = OrderedDict()
			for i,term in enumerate(config_terms):
				hash = ''.join(term)
				if hash in param_sets:
					param_sets[hash][0]+= 1
				else:
					param_sets[hash] = [1,term,config_gibbs[i]]
			
			# order terms by most negative free energy first
			param_sets = OrderedDict(sorted(param_sets.iteritems(),key=lambda x: x[1][2]))
			
			for hash in param_sets:
				if param_sets[hash][1] == []:
					ret+= "%i+"%(param_sets[hash][0])
					continue
					
				ret+= "%i("%(param_sets[hash][0])
				skip = []
				for p in param_sets[hash][1]:
					if p in skip:
						break
					tmp = param_sets[hash][1].count(p)
					if tmp > 0:
						skip.append(p)
						if tmp == 1:
							ret+="%s"%(p)
						else:
							ret+="%s^{%i}"%(p,tmp)
				ret+= ")+"
			
			ret = ret[:-1] # remove last +
			ret+="]\\\\\n"
		
		ret+= "\\end{array}"
		return ret
		
	def set_precision(self,precision=1E-9):
		"""Sets the precision in free ligand concentration for set_probabilities()"""
		self.precision = precision
		
	def set_probabilities(self,totalP,totalL,T):
		"""Set the normalized weights (probabilities) of each configuration at the specified conditions

		Arguments:
			totalP (float): The total concentration of binding (protein) lattices
			totalL (float): The total concentration of ligands
			T (float): The experimental temperature
			
		Returns:
			(float): The free concentration of ligand
		"""

		def _freeL_dev(freeL): # 
			"""Return the deviation between predicted and actual free ligand concentration """
			
			# set the probability of each configuration at the free ligand concentration
			self.weights = [math.exp( (-1.0 * self.gibbs[i]) / ( _R * T ) ) * freeL**self.bound[i] for i in xrange(self.nconfigs) ]
			total = sum(self.weights)
			self.weights = [ self.weights[i] / total for i in xrange(self.nconfigs) ]
			
			# concentration of sites in bound state
			bound = sum( [totalP * self.weights[i] * self.bound[i] for i in xrange(self.nconfigs)] )

			return totalL - (freeL + bound)
		
		# find where the deviation between actual and test free ligand is zero. Use zero free and total ligand as our bracketing guesses
		freeL = scipy.optimize.brentq( _freeL_dev, 0.0, totalL, xtol=self.precision, disp=True )
		
		# make sure before we quit that we set the weights at the correct free ligand conc
		_freeL_dev( freeL )
		
		return freeL
		
	def Q(self,T0,T,concentrations):
		"""Return the enthalpy of the system at each of the specified concentrations """
		
		# set the free energies (and enthalpic energies if necessary) of each configuration
		self.set_energies(T0,T)

		# calculate the enthalpy at each set of conditions
		Q = [0.0]*len(concentrations)
		for i,c in enumerate(concentrations):
			# set the weights (probabilities) of each lattice configuration
			self.set_probabilities(c['Lattice'],c['Ligand'],T)
			
			# enthalpy is sum of all weighted enthalpies of the lattices
			Q[i] = sum( [self.weights[j] * self.enthalpies[j] for j in xrange(self.nconfigs)] )

		return Q

	def set_energies(self,T0,T):
		for i in xrange(self.nconfigs):
			self.gibbs[i], self.enthalpies = 0,0
		raise NotImplementedError("Valid ITC Ising models should implement this!")

class FullAdditive(Ising):
	"""An Ising-type model, in which ligands bind to either a linear or circular lattice.
Coupling can occur to both unoccupied and occupied lattice points."""

	def __init__(self,nsites=3,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding to a site.' )
		self.add_parameter( 'dGa',	'dG',	description='Additional free energy change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dGb',	'dG',	description='Additional free energy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding to a site.' )
		self.add_parameter( 'dHa',	'dH',	description='Additional enthalpy change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dHb',	'dH',	description='Additional enthalpy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic heat capacity change upon binding to a site.' )
		self.add_parameter( 'dCpa',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dCpb',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an occupied site' )

	def set_energies(self,T0,T):
		dG0 = dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 )
		dGa = dG_vant_Hoff( self.params['dGa'], self.params['dHa'], self.params['dCpa'], T, T0 )
		dGb = dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 )
		dH0 = dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 )
		dHa = dH_vant_Hoff( self.params['dHa'], self.params['dCpa'], T, T0 )
		dHb = dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			
			for j in xrange(self.nsites):

				if self.get_site_occupancy(i,j): # is site occupied?
					self.gibbs[i]+=dG0
					self.enthalpies[i]+=dH0
					self.config_params[i].append( 'K_0' )
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						self.gibbs[i]+=dGb
						self.enthalpies[i]+=dGb
						self.config_params[i].append( 'K_b' )

					elif self.circular:
						self.gibbs[i]+=dGa
						self.enthalpies[i]+=dGa
						self.config_params[i].append( 'K_a' )
					
					if self.get_site_occupancy(i,j-1): # is previous neighboring site occupied?
						pass # Note: this avoids double counting, and thus is implemented as in Saroff & Kiefer

					elif self.circular:
						self.gibbs[i]+=dGa
						self.enthalpies[i]+=dGa
						self.config_params[i].append( 'K_a' )
		return

class HalfAdditive(Ising):
	"""An Ising-type model, in which ligands bind to either a linear or circular lattice.
Coupling only occurs between occupied lattice points."""

	def __init__(self,nsites=3,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding to a site.' )
		self.add_parameter( 'dGb',	'dG',	description='Additional free energy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding to a site.' )
		self.add_parameter( 'dHb',	'dH',	description='Additional enthalpy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dCpb',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an occupied site' )

	def set_energies(self,T0,T):
		dG0 = dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 )
		dGb = dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 )
		dH0 = dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 )
		dHb = dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			
			for j in xrange(self.nsites):

				if self.get_site_occupancy(i,j): # is site occupied?
					self.gibbs[i]+=dG0
					self.enthalpies[i]+=dH0
					self.config_params[i].append( 'K_0' )
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						self.gibbs[i]+=dGb
						self.enthalpies[i]+=dGb
						self.config_params[i].append( 'K_b' )
		return

class NonAdditive(Ising):
	"""An Ising-type model, in which ligands bind to either a linear or circular lattice.
Binding thermodynamics depend upon whether zero, one, or both neighboring sites are occupied."""

	def __init__(self,nsites=3,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dGX',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dGY',	'dG',	description='Free energy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dGZ',	'dG',	description='Free energy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dHX',	'dH',	description='Enthalpy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dHY',	'dH',	description='Enthalpy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dHZ',	'dH',	description='Enthalpy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dCpX',	'dCp',	description='Change in heat capacity upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dCpY',	'dCp',	description='Change in heat capacity upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dCpZ',	'dCp',	description='Change in heat capacity upon binding to a site flanked by two occupied' )

	def set_energies(self,T0,T):
		dGX = dG_vant_Hoff( self.params['dGX'], self.params['dHX'], self.params['dCpX'], T, T0 )
		dGY = dG_vant_Hoff( self.params['dGY'], self.params['dHY'], self.params['dCpY'], T, T0 )
		dGZ = dG_vant_Hoff( self.params['dGZ'], self.params['dHZ'], self.params['dCpZ'], T, T0 )
		dHX = dH_vant_Hoff( self.params['dHX'], self.params['dCpX'], T, T0 )
		dHY = dH_vant_Hoff( self.params['dHY'], self.params['dCpY'], T, T0 )
		dHZ = dH_vant_Hoff( self.params['dHZ'], self.params['dCpZ'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			
			for j in xrange(self.nsites):
				if self.get_site_occupancy(i,j): # is site occupied?
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						if self.get_site_occupancy(i,j-1): # is previous neighboring site occupied?
							self.gibbs[i]+= dGZ
							self.enthalpies[i]+= dHZ
							self.config_params[i].append( 'K_Z' )
							
						elif self.circular:
							self.gibbs[i]+= dGY
							self.enthalpies[i]+= dHY
							self.config_params[i].append( 'K_Y' )
							
					elif self.get_site_occupancy(i,j-1):
						self.gibbs[i]+= dGY
						self.enthalpies[i]+= dHY
						self.config_params[i].append( 'K_Y' )
						
					elif self.circular:
						self.gibbs[i]+= dGX
						self.enthalpies[i]+= dHX
						self.config_params[i].append( 'K_X' )
						
		return
