from os		import path
from scipy	import array,dtype,zeros,genfromtxt,std
from thermo	import *
from utilities	import savitzky_golay

class ITCExperiment:
	"""A class that encapsulates a specific ITC experiment.

	Attributes:
		title (string): A descriptive title for the experiment.
		T (float): The experimental temperature (in Kelvin).
		V0 (float): The volume of the calorimeter's cell, in microliters.
		M0 (float): The concentration of the macromolecule used in the experiment.
		L0 (float): The concentration of the ligand used in the experiment.
		dQ_exp (list of floats): The experimental (observed) enthalpies at each injection.
		dQ_fit (list of floats): The predicted enthalpies (if any) at each injection.
		spline (list of floats): The smoothed spline (if any) to the experimental enthalpies.
		I_vol (list of floats): The volumes injected at each titration point.
		M_conc (list of floats): The total macromolecule concentrations in the cell at each titration point.
		L_conc (list of floats): The total ligand concentrations in the cell at each titration point.
		reverse (boolean): Is the experiment performed in reverse? i.e. Macromolecule injected.
		ratios (list of floats): The ratio of syringe component concentration to cell concentration component at each titration point.
		skip (list of ints): Titration points to exclude during fitting.
		dQ_err (list of floats): List of estimated errors in each titration point enthalpy.
		spline_pts (int): Number of points to use in the Savitsky-Golay filter used to estimate errors in experimental enthalpies.
		spline_order (int): Order of the SG filter.
		dil_Q (float): Per-microliter heat of dilution.
		chisq (float): If a fit has been generated, the reduced chi-squared goodness-of-fit value.
	"""

	def __init__(self, title, T, V0, M0, L0, dQ_exp, I_vol, reverse=False, skip=[], dQ_err=[], spline_pts=7, spline_order=1, dil_Q=0.0):
		"""Constructor for the ITCExperiment object.

		Args:
			title (string): A descriptive title for the experiment.
			T (float): The experimental temperature (in Kelvin).
			V0 (float): The volume of the calorimeter's cell, in microliters.
			M0 (float): The concentration of the macromolecule used in the experiment.
			L0 (float): The concentration of the ligand used in the experiment.
			dQ_exp (list of floats): The experimental (observed) enthalpies at each injection.
			I_vol (list of floats): The volumes injected at each titration point.
			reverse (boolean): Is the experiment performed in reverse? i.e. Macromolecule injected.
			skip (list of ints): Titration points to exclude during fitting.
			dQ_err (list of floats): List of estimated errors in each titration point enthalpy.
			spline_pts (int): Number of points to use in the Savitsky-Golay filter used to estimate errors in experimental enthalpies.
			spline_order (int): Order of the SG filter.
			dil_Q (float): Per-microliter heat of dilution.
		"""

		self.npoints = len(dQ_exp)
		assert( self.npoints == len(I_vol) )

		self.title	= title
		self.T		= T
		self.V0		= V0
		self.M0		= M0
		self.L0		= L0
		self.I_vol	= I_vol
		self.reverse= reverse
		self.skip	= skip
		self.dil_Q	= dil_Q

		cM,cL = [0.0]*self.npoints,[0.0]*self.npoints

		# convert raw data (in calories) to joules per mol of injectant
		if self.reverse:
			cL[0] = L0
			self.dQ_exp	= [J_from_cal(dQ_exp[i]/(I_vol[i]*M0)) for i in xrange(self.npoints)]
		else:
			cM[0] = M0
			self.dQ_exp	= [J_from_cal(dQ_exp[i]/(I_vol[i]*L0)) for i in xrange(self.npoints)]

		# calculate the total ligand and macromolecule concentration in the cell at each titration point
		# this uses the dilution formula described in Microcal's data processing in Origin manual
		for i in xrange(0,self.npoints):
			dV = sum(I_vol[0:i])
			if reverse:
				cL[i]	= L0 * ( (1-(dV/(2*V0))) / (1+(dV/(2*V0))) )
				cM_r	= (M0*dV/V0) * (1/(1+(dV/(2*V0)))) # concentration of macromolecule from previous injections
				cM[i]	= (M0*I_vol[i]/V0) + cM_r
			else:
				cM[i]	= M0 * ( (1-(dV/(2*V0))) / (1+(dV/(2*V0))) )
				cL_r	= (L0*dV/V0) * (1/(1+(dV/(2*V0)))) # concentration of ligand from previous injections
				cL[i]	= (L0*I_vol[i]/V0) + cL_r

		# convert the macromolecule and ligand concentrations to C doubles
		self.M_conc	= array(cM,dtype('d'))
		self.L_conc	= array(cL,dtype('d'))

		if self.reverse:
			self.ratios = get_ratios(self.M_conc,self.L_conc)
		else:
			self.ratios = get_ratios(self.L_conc,self.M_conc)

		if dQ_err != []:
			assert len(dQ_err) == self.npoints
			self.dQ_err = dQ_err[:]
			self.spline = None
		else:
			# estimate errors by fitting spline
			counter,tmp = 0,[0]*self.npoints
			for i in xrange(self.npoints):
				if i not in self.skip:
					tmp[i] = counter
					counter += 1

			spl = savitzky_golay([ self.dQ_exp[i] for i in xrange(self.npoints) if i not in self.skip ], spline_pts, spline_order )
			err = std([self.dQ_exp[i] - spl[tmp[i]] for i in xrange(self.npoints) if i not in self.skip])
			self.dQ_err = [err]*self.npoints
			self.spline = [spl[tmp[i]] for i in xrange(self.npoints)]

		self.dQ_fit = [0.0]*self.npoints

	def calc_chisq(self, dQ_fit, writeback=True):
		"""Calculate the goodness-of-fit between the provided data and the experimental data.

		Args:
			dQ_fit (list of floats): The predicted injection enthalpies to compare against the experimental data.
			writeback (boolean): Do we overwrite the object's dQ_fit attribute?

		Returns:
			(float): The goodness of the fit, as a reduced chi-square.
		"""

		if writeback:
			self.dQ_fit = dQ_fit[:]

		self.chisq = 0.0
		for i in xrange(self.npoints):
			if i not in self.skip:
				self.chisq += (self.dQ_exp[i] - dQ_fit[i])**2 / self.dQ_err[i]**2
		self.chisq /= (self.npoints -len(self.skip))
		return self.chisq

	def show_plot(self,hardcopy=False,hardcopydir='.',hardcopyprefix='', hardcopytype='png'):
		"""Generate a plot of the experimental data, and the fit if present.

		Arguments:
			hardcopy (boolean): Display the fit to the screen, or write it to a file?
			hardcopydir (string): The directory to write the hardcopy to.
			hardcopyprefix (string): A prefix to append to the experiment's title when generating the output plot filename.
			hardcopytype (string): The output format for the plot (availability is dependent upon what backends matplotlib was compiled with).

		Returns:
			None

		"""
		try:
			#import matplotlib
			#matplotlib.use('Agg')
			import matplotlib.pyplot as pyplot
		except:
			pyplot = None

		if pyplot == None: return
		if hardcopy: fig = pyplot.figure()

		pyplot.clf()
		pyplot.title(self.title)
		pyplot.ylabel("Joule/mol of injectant")

		if self.reverse:
			pyplot.xlabel("Molar Ratio (M/L)")
		else:
			pyplot.xlabel("Molar Ratio (L/M)")

		tmpx = [ self.ratios[i] for i in xrange(self.npoints) if i not in self.skip ]
		tmpy = [ self.dQ_exp[i] for i in xrange(self.npoints) if i not in self.skip ]
		tmpz = [ self.spline[i] for i in xrange(self.npoints) if i not in self.skip ]
		tmpd = [ self.dQ_err[i] for i in xrange(self.npoints) if i not in self.skip ]

		pyplot.errorbar(tmpx,tmpy,yerr=tmpd,c='#000000',fmt='s')
		if self.spline != None:
			pyplot.plot(tmpx,tmpz,c='g')

		if len(self.dQ_fit) > 1:
			tmpy = [ self.dQ_fit[i] for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.plot(tmpx,tmpy,c='r')

		if len(self.skip) > 0:
			tmpx = [ self.ratios[i] for i in self.skip ]
			tmpy = [ self.dQ_exp[i] for i in self.skip ]
			pyplot.errorbar(tmpx,tmpy,yerr=0,c='g',fmt='s')

		pyplot.draw()
		if hardcopy:
			fig.savefig( path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,self.title,hardcopytype)), bbox_inches='tight')
			pyplot.close(fig)
		else:
			pyplot.show()

	def export_data(self,path,units='J'):
		"""Export the attributes of the experiment to a file.

		Args:
			path (string): The filesystem path to write the output to.
			units (string): The units to use in the export, must be either "J" for Joules (the default), or "cal" for calories.

		Returns:
			None
		"""

		assert units in ('cal','J')

		from datetime import datetime

		h = open(path,'w')
		h.write("#\n#itcsim export of %s\n#\n"%(self.title))
		h.write("#Date	%s\n"%(datetime.ctime(datetime.today())))
		h.write("#T	%.5f\n"%(self.T))
		h.write("#V0	%.5f\n"%(self.V0))
		h.write("#M0	%.5E\n"%(self.M0))
		h.write("#L0	%.5E\n"%(self.L0))
		h.write("#Reverse %s\n"%(str(self.reverse==True)))
		h.write("#dilQ	%.5E\n"%(self.dil_Q))
		h.write("#Skipped: %s\n"%(",".join(map(str,self.skip))))
		h.write("#\n#	Ivol	M	L	Ratio	dQ_exp	dQ_fit	spline	ddQ_exp	skipped\n")

		if self.spline == None:
			spline = [0 for i in xrange(self.npoints)]
		else:
			spline = self.spline[:]

		def _format(f):
			if units == 'cal':
				return J_to_cal(f)
			return f

		for i in xrange(self.npoints):
			h.write("%i	%.5f	%.5E	%.5E	%.5f	%.5E	%.5E	%.5E	%.5E	%s\n"%(
				i,
				self.I_vol[i],
				self.M_conc[i],
				self.L_conc[i],
				self.ratios[i],
				_format(self.dQ_exp[i]),
				_format(self.dQ_fit[i]),
				_format(spline[i]),
				_format(self.dQ_err[i]),
				(i in self.skip)
			))

		h.close()
