"""Extension of Ising classes for fitting mass spectrometry data.


"""

import os
import sys
import scipy
import uuid

from .itc_experiment import ITCExperimentBase
from .model_ising import Ising


_MATPLOTLIB_BACKEND = None #None for default
_FILE_COMMENT = ["#"]
_FILE_KEYPAIR = ["@"]

class MSExperiment(ITCExperimentBase):
	"""A class that represents a specific mass spec experiment, i.e. weights of different stoichiometry populations

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

	def __init__(self, path, T=298.15, sigma=0.05, title=None, lattice_name="Lattice", ligand_name="Ligand"):
		"""Constructor for the MSExperiment object.

		Arguments
		---------
		path : string
			A path to a file containing experimental populations in the format "[Lattice] [Ligand] [population state 0] [population state 1]..."
		T : float
			An uncertainty in the measured population state abundances (default is 5%).
		sigma : float
			An uncertainty in the measured population state abundances (default is 5%). If None, read uncertainties from the experimental file.
		title : string
			A descriptive title for the experiment (default is the trimmed filename).
		"""
		
		assert os.path.isfile(path)

		# relevant file parameters. TODO: replace using getattr with defaults?
		keypairs = {
			"Lattice":lattice_name,
			"Ligand":ligand_name,
			"Temperature":T,
			"Error":sigma,
			"Title":title,
		}

		# dummy vars to execute the ITCExperimentBase base class constructor.
		V0, injections, dQ = 1.0, [1.0], [1.0]
		Cell, Syringe = {keypairs['Lattice']:1.0},{keypairs['Ligand']:1.0}
		ITCExperimentBase.__init__(self, keypairs['Temperature'], V0, injections, dQ, Cell, Syringe)

		# reset relevant attributes
		self.path = path
		self.Concentrations = []
		self.npops, self.npoints = None, 0

		data = []
		with open(path) as fh:
			for i, line in enumerate(fh):

				# ignore empty and comment lines
				if line.strip() == "":
					continue
				elif line[0] in _FILE_COMMENT: # comment
					if line[1] in _FILE_KEYPAIR: # convert "#@" to "@" to play nice with some plotting scripts
						line = line[1:]
					else:
						continue

				if line[0] in _FILE_KEYPAIR: # set a variable
					tmp = line.split()
					try:
						assert len(tmp) > 2
					except AssertionError:
						raise Exception("Found malformed experiment parameter specification at line %i, should be of the form \"@param = value\""%(i+1))

					try:
						assert tmp[0][1:] in keypairs
					except AssertionError:
						raise Exception("Found unknown parameter specification \"%s\" at line %i."%(tmp[0][1:],i+1))
					
					try:
						if type(keypairs[tmp[0][1:]]) is str or keypairs[tmp[0][1:]] is None: # title
							keypairs[tmp[0][1:]] = " ".join(tmp[2:])
						else:
							keypairs[tmp[0][1:]] = float(tmp[2])
					except ValueError:
						raise Exception("Found malformed parameter value (could not convert \"%s\" to float) at line %i"%(tmp[2],i+1))
					continue

				arr = line.split()
				if len(arr) < 3:
					raise Exception("Found malformed titration point (less than three columns) at line %i"%(i+1))		
				elif self.npops == None:
					self.npops = len(arr) -2
				elif self.npops != len(arr) -2:
					raise Exception("Found malformed titration point (inconsistent number of columns) at line %i"%(i+1))
				
				try:
					tmp = list(map(float,arr))
				except ValueError:
					raise Exception("Found malformed titration point (could not convert column values to floats) at line %i"%(i+1))

				data.extend(tmp[2:]) # strip the lattice/ligand concentration columns

				self.Concentrations.append({})
				self.Concentrations[-1][keypairs['Lattice']] = tmp[0]
				self.Concentrations[-1][keypairs['Ligand']] = tmp[1]

				self.npoints += self.npops

		if keypairs['Title'] is None:
			self.title = os.path.splitext(os.path.basename(self.path))[0]
		else:
			self.title = keypairs['Title']
		
		self.lattice_name, self.ligand_name = keypairs['Lattice'], keypairs['Ligand']

		if keypairs['Error'] is None:
			try:
				assert self.npoints % 2 == 0
				assert min(data[int(self.npoints/2):]) > 0.0
			except AssertionError:
				raise Exception("If sigma == None, must provide nonzero experimental uncertainties in the experimental file.")
				
			self.npoints = int(self.npoints/2) # first half of points matrix are experimental values, second half are uncertainties
			self.Concentrations = self.Concentrations[:self.npoints] # discard the duplicated concentration columns for sigmas
			self.PopIntens = scipy.array(data[:self.npoints]).reshape((int(self.npoints/self.npops),self.npops))
			self.PopSigmas = scipy.array(data[self.npoints:]).reshape((int(self.npoints/self.npops),self.npops))
		else:
			self.sigma = keypairs['Error']
			self.PopIntens = scipy.array(data).reshape((int(self.npoints/self.npops),self.npops))
			self.PopSigmas = scipy.full(self.PopIntens.shape,self.sigma)

		# normalize intensities to 1, accordingly scale their sigmas
		totals = self.PopIntens.sum(axis=1,keepdims=True)
		self.PopIntens /= totals
		self.PopSigmas /= totals

		# precompute variance (s**2) from sigmas
		self.PopSigmas = scipy.square(self.PopSigmas)
		self.sigma = scipy.sqrt(scipy.mean(self.PopSigmas))

		self.PopFits = scipy.zeros(self.PopIntens.shape)

		self.chisq = None
		self.initialized = True

	def __str__(self):
		"""Stringify the experiment to be suitable for display to the user
		
		Parameters
		----------
		None
		
		Returns
		-------
		string
			String containing the description and other experimental info.
		"""
		
		ret = "Title: %s\n"%(self.title)
		if (self.chisq != None):
			ret+= "Current chisq: %f\n"%(self.chisq)
		ret+= "Temperature: %0.2f K\n"%(self.T)
		ret+= "Titration points: %i\n"%(self.npoints/self.npops)
		ret+= "Stoichiometries: %i\n"%(self.npops)
		ret+= "Components:\n"
		ret+= "\tLattice: \"%s\"\n"%(self.lattice_name)
		ret+= "\tLigand: \"%s\"\n"%(self.ligand_name)
		return ret
	
	def make_plot(self, hardcopy=False, hardcopydir='.', hardcopyprefix='', hardcopytype='png'):
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

		if _MATPLOTLIB_BACKEND != None:
			import matplotlib
			matplotlib.use(_MATPLOTLIB_BACKEND)
		import matplotlib.pyplot as pyplot

		if hardcopy: fig = pyplot.figure()

		pyplot.clf()
		pyplot.title(self.title)
		pyplot.ylabel("Experimental - Fit Abundance")
		pyplot.xlabel("Mass Populations")
		
		ax1 = pyplot.subplot(3,1,1)
		ax2 = pyplot.subplot(3,1,2)
		ax3 = pyplot.subplot(3,1,3)

		x_points = int(self.npoints/self.npops)

		xax_positions,xax_labels = [],[]
		width,space,left = 0.25,0.5,0.0
		for i in range(self.npops):
			bars1 = ax1.bar( [left + (j*width) for j in range(x_points)], [self.PopIntens[j][i] for j in range(x_points)], width=width, edgecolor='r' )
			bars2 = ax2.bar( [left + (j*width) for j in range(x_points)], [self.PopFits[j][i] for j in range(x_points)], width=width, edgecolor='r' )
			bars3 = ax3.bar( [left + (j*width) for j in range(x_points)], [self.PopIntens[j][i] - self.PopFits[j][i] for j in range(x_points)], width=width, edgecolor='r' )
			xax_positions.append( left + (x_points*width)/2.0 )
			xax_labels.append("%i"%i)
			left += (x_points*width)+space
	
			for j in range(x_points):
				color = 1.0 - (float(j) / x_points)
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
			fig.savefig( os.path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,self.title,hardcopytype)), bbox_inches='tight')
			pyplot.close(fig)
		else:
			pyplot.show()

	def make_population_plot(self, dataset='fit', hardcopy=False, hardcopydir='.', hardcopyprefix='', hardcopytype='png'):
		assert dataset in ("fit", "experimental", "residuals")

		if not self.initialized and (dataset=="fit" or dataset=="residuals"):
			raise AssertionError("Fit data has not yet been generated for this experiment")

		if _MATPLOTLIB_BACKEND != None:
			import matplotlib
			matplotlib.use(_MATPLOTLIB_BACKEND)
		
		import matplotlib.pyplot as pyplot
		from mpl_toolkits.mplot3d import axes3d, Axes3D
		
		pyplot.clf()
		fig = pyplot.figure()
		ax = Axes3D(fig)
#		ax = fig.add_subplot(1, 1, 1, projection='3d') matplotlib version issue

		xs = range(int(self.npoints / self.npops)) # number of concentrations
		for i in range(self.npops): # iterating over each configuration 

			ys_fit = scipy.array([pop[i] for pop in self.PopFits])
			ys_exp = scipy.array([pop[i] for pop in self.PopIntens])
			
			if dataset == "fit":
				ys = ys_fit
			elif dataset == "experimental":
				ys = ys_exp
			else:
				ys = ys_exp - ys_fit

			# z-order doesn't actually do anything afaik, it's just here to show I tried
			ax.bar(xs, ys, zs=i, zdir='y', zorder=(self.npops-i), color=pyplot.cm.jet(1.0 * i / self.npops), alpha=0.8)

		ax.set_title("%s:(%s)"%(self.title, dataset))
		ax.set_xlabel('Titration Point')
		ax.set_ylabel('Stoichiometry')
		ax.set_zlabel('Abundance')
		ax.set_zlim([0,1])
		pyplot.draw()

		if hardcopy:
			fig.savefig( os.path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,self.title,hardcopytype)), bbox_inches='tight')
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
		for i in range(int(self.npoints/self.npops)):
			fh.write("%0.3E\t%0.3E\t%s\n" % (
				self.Concentrations[i][self.lattice_name],
				self.Concentrations[i][self.ligand_name],
				"\t".join(["%f"%f for f in self.PopIntens[i]])))

		fh.write("# Fitted data\n")
		for i in range(int(self.npoints/self.npops)):
			fh.write("%0.3E\t%0.3E\t%s\n" % (
				self.Concentrations[i][self.lattice_name],
				self.Concentrations[i][self.ligand_name],
				"\t".join(["%f"%f for f in self.PopFits[i]])))

		fh.close()
			
	def get_chisq(self, pops, writeback=False):
		"""Calculate the chi-square goodness-of-fit between the experimental population abundances and the fitted ones.

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
			Variances (sigma**2) must have been precomputed as self.PopSigmas (perhaps should be self.PopVariances?)
		"""

		assert self.PopIntens.shape == pops.shape

		self.chisq = scipy.sum(scipy.square(self.PopIntens - pops) / self.PopSigmas) / self.PopIntens.size

		if writeback:
			self.PopFits = pops		

		return self.chisq

class MSExperimentSynthetic(MSExperiment):
	"""A MSExperiment-derived class that can be used to simulate a mass spec experiment"""
		
	def __init__(self, lattice_concs, ligand_concs, T=298.15, noise=0.02, title=None, lattice_name="Lattice", ligand_name="Ligand"):
		"""Constructor for the MSExperimentSynthetic object.

		Arguments
		---------
		lattice_concentrations : list of floats
			Concentrations of the lattice component (in M)
		lattice_concentrations : list of floats
			Concentrations of the ligand component (in M)
		T : float
			The temperature of the experiment to simulate.
		noise : string
			The magnitude of noise to add to the simulated data.
		title :
			The title to assign to the experiment
		"""

		assert len(lattice_concs) == len(ligand_concs)

		self.sigma = noise
		if title is not None:
			self.title = title
		else:
			self.title = uuid.uuid4()

		self.lattice_name, self.ligand_name = lattice_name, ligand_name

		# patches to play nice with the ITCExperimentBase base class constructor. Do we even need this?
		V0,injections,dQ = 1.0,[1.0],[1.0]
		Cell,Syringe = {self.lattice_name:1.0},{self.ligand_name:1.0}
		self.path = None

		ITCExperimentBase.__init__(self, T, V0, injections, dQ, Cell, Syringe, title=self.title)

		# reset key attributes now
		self.Concentrations = [{self.lattice_name:lattice_concs[i],self.ligand_name:ligand_concs[i]} for i in range(len(lattice_concs))]

		self.initialized = False
		self.chisq = None

	def __str__(self):
		if not self.initialized:
			return "Synthetic experiment has not yet been generated - attach to a simulator and run to populate"
		return MSExperiment.__str__(self)

	def get_chisq(self, pops, writeback=False):
		"""Extends the base MSExperiment class to overwrite the PopFits attribute on the first call"""

		if not self.initialized:
			self.PopIntens = scipy.absolute(scipy.random.normal(pops,self.sigma))
			self.PopSigmas = scipy.full(self.PopIntens.shape,self.sigma**2)
			self.PopFits = scipy.zeros(self.PopIntens.shape)
			self.npoints, self.npops = self.PopIntens.shape
			self.npoints = self.npoints*self.npops
			self.initialized = True

		return MSExperiment.get_chisq(self, pops, writeback)


class MSModel(Ising):
	"""
	Class used to convert an existing Ising-based binding model to one that returns lattice+ligand stoichiometry abundances suitable for fitting mass spec data.
	
	Note #1: Inherits all attributes of the provided class instance, serves as a passthru using the provided class' set_energies().
	Note #2: If the provided class instance has overwritten any of the Ising class methods (other than set_energies), you're best off not using this convertor class."""

	def __init__(self,model):
		
		def _recursor(base, base_list):
			if base.__bases__ == ():
				return []
			else:
				for b in base.__bases__:
					base_list.append( b.__name__ )
					_recursor( b, base_list)
				return base_list

		# ensure that an Ising-based model has been in fact passed
		assert "Ising" in _recursor(model.__class__, [])

		# copy all of the (initialized) parent model attributes
		self.__dict__ = model.__dict__.copy()
		self.model = model

		self.__doc__ = "A model adapted from base type \"%s.%s\" for fitting mass spectrometric population data.\n\nOriginal docstring:\n%s"%(self.model.__module__,self.model.__class__.__name__,self.model.__doc__) 

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
		
		ret = scipy.zeros((len(concentrations),self.model.nsites+1))
		for i,c in enumerate(concentrations):
			
			# set the probabilities (weights) for all configurations
			self.model.set_probabilities(c[self.lattice_name],c[self.ligand_name],T)
			
			# for each stoichiometry, add up the configuration weights
			for j in range(self.model.nconfigs):
				ret[i][self.model.bound[j]] += self.model.weights[j]

		return ret