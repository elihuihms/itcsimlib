import os
import sys
import numpy as np

from itc_experiment import ITCExperimentBase
from model_ising import Ising

_MATPLOTLIB_BACKEND = None #None for default

class MSExperiment(ITCExperimentBase):
	"""A class that encapsulates a specific mass spec epxeriment, i.e. weights of different stoichiometry populations

	Attributes
	----------
	chisq : float
		If a fit has been generated, the reduced chi-squared goodness-of-fit value.
	T : float
		The experimental temperature (in Kelvin).
	Lattice : string
		The name of the lattice component (if it exists, read from the experiment file)
	Ligand : string
		The name of the ligand component (if it exists, read from the experiment file)
	Concentrations : list of dicts
		The concentrations of the components at each titration point
	PopIntens : ndarray
		The normalized intensities of each lattice+ligand stoichiometry at each set of component concentrations
	PopSigmas : ndarray
		The uncertainty in the measured intensity of each lattice+ligand stoichiometry at each set of component concentrations
	PopFits	: ndarray
		The fitted intensities (if calculated) of each lattice+ligand stoichiometry at each set of component concentrations

	Notes
	-----
		This class abuses several of the usual ITCExperimentBase attributes in a rather dirty fashion - it won't as a simple replacement for many itcsimlib routines, but basic fitting should be OK
	"""
		
	def __init__(self, path, T=298.15, sigma=0.05, title=None):
		"""Constructor for the MSExperiment object.

		Arguments
		---------
		path : string
			A path to a file containing experimental populations in the format "[Lattice] [Ligand] [population state 0] [population state 1]..."
		T : float
			An uncertainty in the measured population state abundances (default is 5%).
		sigma : float
			An uncertainty in the measured population state abundances (default is 5%).
		title : string
			A descriptive title for the experiment (default is the trimmed filename).
		"""
		
		assert os.path.isfile(path)
		assert sigma > 0.0

		self.sigma = sigma
		
		self.path = path
		if title:
			self.title = title
		else:
			self.title = os.path.splitext(os.path.basename(self.path))[0]

		# patches to play nice with the ITCExperimentBase base class constructor
		V0,injections,dQ = 1.0,[1.0],[1.0]
		Cell,Syringe = {"Lattice":1.0},{"Ligand":1.0}

		ITCExperimentBase.__init__(self, T, V0, injections, dQ, Cell, Syringe, title=self.title)

		# reset key attributes now
		self.Concentrations = []
		self.npops, self.npoints = None, 0
		self.Lattice, self.Ligand = "Lattice", "Ligand"

		data = []
		with open(path) as fh:
			for i, line in enumerate(fh):

				# ignore empty and comment lines
				if line.strip() == "":
					continue
				elif line[0] in ["#"]: # comment
					continue
				elif line[0] in ["@"]: # set a variable
					tmp = line.split()

					# hacky experimental file parameter parser. TODO: use ConfigParser in Python 3 to handle file header string
					try:
						if tmp[0] == "@Temperature":
							try:
								self.T = float(tmp[2])
							except ValueError:
								raise Exception("Found malformed Temperature argument in experiment header (could not convert \"%s\" to a number) at line %i"%(tmp[2],i+1))
						elif tmp[0] == "@Error":
							try:
								self.sigma = float(tmp[2])
							except ValueError:
								raise Exception("Found malformed Error argument in experiment header (could not convert \"%s\" to a number) at line %i"%(tmp[2],i+1))
						elif tmp[0] == "@Title":
							self.title = " ".join(tmp[2:])
						elif tmp[0] == "@Lattice":
							self.Lattice = tmp[2]
						elif tmp[0] == "@Ligand":
							self.Ligand = tmp[2]
					except IndexError:
						raise Exception("Found malformed experiment parameter specification at line %i, should be of the format \"@param = value\""%(i+1))

					continue

				# basic error checking
				arr = line.split()
				if len(arr) < 3:
					raise Exception("Found malformed titration point (less than three columns) at line %i"%(i+1))		
				elif self.npops == None:
					self.npops = len(arr) -2
				elif self.npops != len(arr) -2:
					raise Exception("Found malformed titration point (inconsistent number of columns) at line %i"%(i+1))
				
				try:
					tmp = map(float,arr)
				except ValueError:
					raise Exception("Found malformed titration point (could not convert column values to floats) at line %i"%(i+1))

				data.extend([f/sum(tmp[2:]) for f in tmp[2:]]) # append the normalized abundances at each titration point
				self.Concentrations.append({})
				self.Concentrations[-1]["Lattice"]	= tmp[0] / 1.0E6
				self.Concentrations[-1]["Ligand"]	= tmp[1] / 1.0E6

				self.npoints += 1

		self.PopIntens	= np.array(data).reshape((self.npoints,self.npops))
		self.PopSigmas	= np.full(self.PopIntens.shape,self.sigma**2)
		self.PopFits	= np.zeros(self.PopIntens.shape)

		self.chisq = None
	
	def make_plot(self,hardcopy=False,hardcopydir='.',hardcopyprefix='', hardcopytype='png'):
		"""Generate a stacked plot of the experimental populations, the fitted populations, and the residuals between the two.

		Arguments
		----------
		hardcopy : boolean
			Display the fit to the screen, or write it to a file?
		hardcopydir : string
			The directory to write the hardcopy to.
		hardcopyprefix : string
			A prefix to append to the experiment's title when generating the output plot filename.
		hardcopytype : string
			The output format for the plot (availability is dependent upon what backends matplotlib was compiled with).

		Returns
		-------
		None
		"""

		try:
			if _MATPLOTLIB_BACKEND != None:
				import matplotlib
				matplotlib.use(_MATPLOTLIB_BACKEND)
			import matplotlib.pyplot as pyplot
		except:
			pyplot = None

		if pyplot == None: return
		if hardcopy: fig = pyplot.figure()

		pyplot.clf()
		pyplot.title(self.title)
		pyplot.ylabel("Experimental - Fit Abundance")
		pyplot.xlabel("Mass Populations")
		
		ax1 = pyplot.subplot(3,1,1)
		ax2 = pyplot.subplot(3,1,2)
		ax3 = pyplot.subplot(3,1,3)

		xax_positions,xax_labels = [],[]
		width,space,left = 0.25,0.5,0.0
		for i in xrange(self.npops):
			bars1 = ax1.bar( [left + (j*width) for j in xrange(self.npoints)], [self.PopIntens[j][i] for j in xrange(self.npoints)], width=width, edgecolor='r' )
			bars2 = ax2.bar( [left + (j*width) for j in xrange(self.npoints)], [self.PopFits[j][i] for j in xrange(self.npoints)], width=width, edgecolor='r' )
			bars3 = ax3.bar( [left + (j*width) for j in xrange(self.npoints)], [self.PopIntens[j][i] - self.PopFits[j][i] for j in xrange(self.npoints)], width=width, edgecolor='r' )
			xax_positions.append( left + (self.npoints*width)/2.0 )
			xax_labels.append("%i"%i)
			left += (j*width)+space
	
			for j in xrange(self.npoints):
				color = 1.0 - (float(j) / self.npoints)
				bars1[j].set_color( (color, color, 1) )	
				bars2[j].set_color( (color, 1, color) )	
				bars3[j].set_color( (1, color, color) )

				bars1[j].set_edgecolor( (0,0,0) )	
				bars2[j].set_edgecolor( (0,0,0) )	
				bars3[j].set_edgecolor( (0,0,0) )	
		
		ax1.set_ylabel("Experimental")
		ax1.set_xticks( xax_positions )
		ax1.set_xticklabels( xax_labels)
		ax2.set_ylabel("Fit")
		ax2.set_xticks( xax_positions )
		ax2.set_xticklabels( xax_labels)
		ax3.set_ylabel("Residuals")
		ax3.set_xticks( xax_positions )
		ax3.set_xticklabels( xax_labels)
			
		pyplot.draw()

		if hardcopy:
			tmp = os.path.splitext(os.path.basename(self.path))[0]
			fig.savefig( os.path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,tmp,hardcopytype)), bbox_inches='tight')
			pyplot.close(fig)
		else:
			pyplot.show()
			
	def export_to_file(self, path):
		"""Export the experimental data and any fit to the specified file path.

		Arguments
		---------
		path : string
			The filesystem path to write the output to.

		Returns
		-------
		None			
		"""
		fh = open(path, 'w')
		
		fh.write("# %s"%self.title)
		fh.write("# Experimental data\n")
		for i in xrange(self.npoints):
			fh.write("%f\t%f\t%s\n" % (
				self.Concentrations[i]['Lattice']*1E6,
				self.Concentrations[i]['Ligand']*1E6,
				"\t".join(["%f"%f for f in self.PopIntens[i]])))

		fh.write("# Fitted data\n")
		for i in xrange(self.npoints):
			fh.write("%f\t%f\t%s\n" % (
				self.Concentrations[i]['Lattice']*1E6,
				self.Concentrations[i]['Ligand']*1E6,
				"\t".join(["%f"%f for f in self.PopFits[i]])))

		fh.close()
			
	def get_chisq(self, pops, writeback=False):
		"""Calculate the goodness-of-fit between the experimental population abundances and the fitted ones.

		Arguments
		---------
		pops : ndarray
			The normalized abundances of each lattice+ligand stoichiometries at each of the provided component concentrations.
		writeback : boolean
			Update the experiment's simulated heat attribute (dQ_fit) with the provided Qs?
			
		Returns
		-------
		float
			The goodness of the fit, as a reduced chi-square.

		Notes
		-----
			This method takes advantage of the fact that ITCSim doesn't inspect the data that is returned by the model, and instead lets the associated experiment handle the goodness-of-fit.
			The only caveat of course is that model must return the same number of lattice+ligand stoichiometries as are present in the experimental results.
		"""

		assert self.PopIntens.shape == pops.shape

		self.chisq = np.sum(np.square(self.PopIntens - pops) / self.PopSigmas) / self.PopIntens.size

		if writeback:
			self.PopFits = pops		

		return self.chisq

class MSModel(Ising):
	"""
	Class used to convert an existing Ising-based binding model to one that returns lattice+ligand stoichiometry abundances suitable for fitting mass spec data.
	
	Note #1: Inherits all attributes of the provided class instance, serves as a passthru using the provided class' set_energies().
	Note #2: If the provided class instance has overwritten any of the Ising class methods (other than set_energies), you're best off not using this convertor class."""

	def __init__(self,model):
		# ensure that an Ising-based model has been in fact passed
		assert "Ising" in [b.__name__ for b in model.__class__.__bases__]

		# copy all of the (initialized) parent model attributes
		self.__dict__ = model.__dict__.copy()
		self.model = model

	def set_energies(self,T0,T):
		"""Update the parent model parameters with whatever we currently have, and set the parent model config energies"""
		self.units = self.model.units
		self.model.params = self.params
		self.model.set_energies(T0,T)

	def Q(self,T0,T,concentrations):
		"""Return a 2D numpy array consisting of the base model's relative stoichiometries at each of the provided component concentrations.
		
		Arguments
		---------
		T0 : float
			The reference temperature of the simulation.
		T : float
			The temperature of the experiment to simulate.
		concentrations : list of dicts
			The concentrations of each component at each titration point.
		
		Returns
		-------
		ndarray
			The normalized abundances of each lattice+ligand stoichiometries at each of the provided component concentrations.
		"""

		# set the energies of this model's configs from the base model
		self.set_energies(T0,T)

		ret = np.zeros((len(concentrations),self.model.nsites+1))
		for i,c in enumerate(concentrations):
			
			# set the probabilities (weights) for all configurations
			self.model.set_probabilities(c['Lattice'],c['Ligand'],T)
			
			# for each stoichiometry, add up the configuration weights
			for j in xrange(self.model.nconfigs):
				ret[i][self.model.bound[j]] += self.model.weights[j]

		return ret