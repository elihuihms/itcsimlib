import os
try:
	import matplotlib.pyplot as pyplot
except:
	pyplot = None

from numpy				import array,dtype
from scipy				import genfromtxt
from multiprocessing	import cpu_count

import itc_lib
from thermo_functions	import *

class ITCExperiment:
	"""
	Encapsulates all the necessary data for a specific ITC experiment
	"""

	def __init__(
		self,
		title,
		T,
		V0,
		M0,
		L0,
		dQ_exp,
		I_vol,
		M_conc,
		L_conc,
		reverse,
		skip=[],
		dG=[],
		dH=[],
		dil_Q=0.0 ):

		assert( len(dQ_exp) == len(I_vol) )
		assert( len(dQ_exp) == len(M_conc) )
		assert( len(dQ_exp) == len(L_conc) )
		assert( len(dG) == len(dH) )

		self.title	= title
		self.T		= T
		self.V0		= V0
		self.P0		= M0
		self.L0		= L0
		self.I_vol	= I_vol
		self.M_conc	= M_conc
		self.L_conc	= L_conc
		self.reverse= reverse
		self.skip	= skip
		self.dil_Q	= dil_Q

		self.dQ_exp = dQ_exp
		self.dQ_fit = []
		self.dG		= dG
		self.dH		= dH
		self.fig = None
		return

	def show_plot(self,hardcopy=False,hardcopydir='.'):
		"""
		Generate a plot of the experimental data, and a fit if available
		"""

		if pyplot == None:
			return

		if self.fig == None and hardcopy:
			self.fig = pyplot.figure()

		pyplot.clf()
		tmp = [ [],[] ]
		for (i,dH) in enumerate(self.dQ_exp):
			if i in self.skip:
				tmp[0].append( self.L_conc[i]/self.M_conc[i] )
				tmp[1].append( dH )

		pyplot.title(self.title)
		pyplot.scatter(get_ratios(self.L_conc,self.M_conc),self.dQ_exp,c='#DDDDDD')
		if len(self.skip) > 0:
			pyplot.scatter(tmp[0],tmp[1],c='g')
		if len(self.dQ_fit) > 1:
			pyplot.plot(get_ratios(self.L_conc,self.M_conc)[1:],self.dQ_fit[1:],c='r')

		#pyplot.ion()
		pyplot.draw()
		if hardcopy:
			self.fig.savefig( os.path.join(hardcopydir,"%s.png"%(self.title)), bbox_inches='tight')
			pyplot.close(self.fig)
		else:
			pyplot.show()

		return

class ITCSim:
	"""
	Contains all the necessary data and methods to collect experiments
	"""

	def __init__(self,T_ref=298.15,model='sk',size=11,cyclic=False,threads=None):
		self.T_ref = T_ref # reference temperature
		self.model = model
		self.size = size
		self.cyclic = cyclic

		if(self.model not in ['simple','sk','ian','hybrid','jump','jump2','jumphybrid']):
			print "Unrecognized model."
			return
		self.lib = "itc_sim_%s.o" % (self.model)

		self.experiments = {} # temperature-keyed dict containing datasets
		return

	def set_params(self, dG, dH, dCp ):
		"""
		Set the various initial parameters for the specified model
		"""

		# note, K & Q values are at ref temp
		self.dG		=	list(dG)
		self.dH		=	list(dH)
		self.dCp	=	list(dCp)

		for T in self.experiments.keys(): # temperatures
			for E in self.experiments[T]:
				E.dG = [dG_vant_Hoff( self.dG[i], self.dH[i], self.dCp[i], T, self.T_ref ) for i in xrange(len(self.dG))],
				E.dH = [dH_vant_Hoff( self.dH[i], self.dCp[i], T, self.T_ref ) for i in xrange(len(self.dG))],

		return

	def add_experiment( self, title, T, V0, M0, L0, DQ, I_vol, reverse=False, skip=[], dil_Q=0.0 ):
		"""
		Add ITC data, along with the necessary experimental and instrument parameters necessary to generate fits
		Calculates the active concentration of protein and ligand at each point
		"""

		if T not in self.experiments.keys():
			self.experiments[T] = []
			self.experiments[T] = []

		cM,cL = [0.0]*len(DQ),[0.0]*len(DQ)

		if reverse:
			cL[0] = L0
		else:
			cM[0] = M0

		for i in range(0,len(DQ)):
			dV = sum(I_vol[0:i])
			if reverse:
				DQ[i] /= (I_vol[i]*M0) # normalize DQ per mol of injectant
				cL[i]	= L0 * ( (1-(dV/(2*V0))) / (1+(dV/(2*V0))) )
				cM_r	= (M0*dV/V0) * (1/(1+(dV/(2*V0)))) # concentration of macromolecule from previous injections
				cM[i]	= (M0*I_vol[i]/V0) + cM_r
			else:
				DQ[i] /= (I_vol[i]*L0)
				cM[i]	= M0 * ( (1-(dV/(2*V0))) / (1+(dV/(2*V0))) )
				cL_r	= (L0*dV/V0) * (1/(1+(dV/(2*V0)))) # concentration of ligand from previous injections
				cL[i]	= (L0*I_vol[i]/V0) + cL_r

				# Ian's calculation of current ligand concentration. As footnoted in his thesis, this is INCORRECT!
				"""
				# % Concentration of injectant in cell assuming RID
				if i>0:
					cL_r = ( cL[i-1] * (V0-I_vol[i]) + (L0*I_vol[i]) ) / V0;
				else:
					cL_r = ( 0.00000 * (V0-I_vol[i]) + (L0*I_vol[i]) ) / V0;
				# % This value is reduced due to some displacement of injected material
				cL[i]	= cL_r*(1-I_vol[i]/2.0/V0)
				"""

		self.experiments[T].append(
			ITCExperiment(
				title	= title,
				T		= T,
				V0		= V0,
				M0		= M0,
				L0		= L0,
				dQ_exp	= array(DQ,dtype('d')),
				I_vol	= array(I_vol,dtype('d')),
				M_conc	= array(cM,dtype('d')),
				L_conc	= array(cL,dtype('d')),
				reverse	= reverse,
				skip	= skip,
				dG		= [dG_vant_Hoff( self.dG[i], self.dH[i], self.dCp[i], T, self.T_ref ) for i in xrange(len(self.dG))],
				dH		= [dH_vant_Hoff( self.dH[i], self.dCp[i], T, self.T_ref ) for i in xrange(len(self.dG))],
				dil_Q	= dil_Q
			)
		)

		return

	def load_file( self, path, T, V0, M0, L0, reverse=False, skip=[], dil_Q=0.0, title=None ):
		"""
		Read experimental ITC data from a file containing two columns: the observed deltaH and the injection volumes
		"""

		(DQ,I_vol) = genfromtxt( path, unpack=True, usecols=(0,1) )
		if title==None:
			title = os.path.splitext(os.path.basename(path))[0]
		self.add_experiment( title, T, V0, M0, L0, DQ, I_vol, reverse, skip, dil_Q )
		return

	def make_plots(self,indices=None,hardcopy=False,hardcopydir='.'):
		"""
		Generate plots for all experimental datasets
		"""

		for T in self.experiments.keys(): # temperatures
			for (i,E) in enumerate(self.experiments[T]):
				if(indices==None) or (i in indices):
					E.show_plot(hardcopy,hardcopydir)
		return

	def update_params( self ):

		for T in self.experiments.keys():
			for E in self.experiments[T]:
				E.dG = dG_vant_Hoff( self.dG[i], self.dH[i], self.dCp[i], T, self.T_ref )
				E.dH = dH_vant_Hoff( self.dH[i], self.dCp[i], T, self.T_ref )

		return

	def make_fits( self,  ):
		"""
		Using the current model parameters, generate fits for all experimental datasets
		"""

		n = len(self.dG)
		params = [0.0]*(2*n)

		ret = {}
		for T in self.experiments.keys(): # temperatures
			for E in self.experiments[T]:

				for i in range(n):
					params[i]	= E.dG[i]
					params[i+n]	= E.dH[i]

				ret[E.title] = itc_lib.sse(
					array(params,dtype('d')),
					self,
					E,
					writeback=True)

		return ret
