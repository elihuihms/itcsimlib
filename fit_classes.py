from multiprocessing	import Pool,cpu_count
from scipy				import optimize

import fit_functions
import opt_functions

class ITCFit:

	def __init__(self,threads=None,xtol=1,ftol=1):
		self.threads = threads # number of processing threads to use

		self.xtol = xtol # convergence tolerance
		self.ftol = ftol

		if(threads!=None):
			self.threads = threads
		else:
			self.threads = cpu_count()

		return

	def fit_H(self, sim):
		"""
		Find optimal dH values for all experimental datasets, without touching model dG values
		"""

		if len(sim.experiments) > 1:
			print "Error:   Cannot fit experiments at more than one temperature with just K"
			return
		if sim.T_ref != sim.experiments.keys()[0] :
			print "Warning: Experimental data is at different temperature (%.2fK) than T_ref (%.2fK)!" % (sim.experiments.keys()[0],sim.T_ref)

		x0 = sim.dH[:]
		ret = optimize.fmin(
			fit_functions.fit_H,
			x0=x0,
			args=(self,sim),
			xtol=self.xtol,
			ftol=self.ftol,
			disp=True)

		print "Fitting of dH concluded:"
		print " dH0 parameters:"
		for (i,f) in enumerate(ret):
			print "  %i : %.5E" % (i,f)

		return ret

	def fit_GH(self, sim, writeback=False):
		"""
		Find optimal dG and dH values for all experimental datasets
		"""

		if len(sim.experiments) > 1:
			print "Error:   Cannot fit experiments at more than one temperature with just K & Q"
			return
		if sim.T_ref != sim.experiments.keys()[0] :
			print "Warning: Experimental data is at different temperature (%.2fK) than T_ref (%.2fK)!" % (sim.experiments.keys()[0],sim.T_ref)

		x0 = list(sim.dG) + list(sim.dH)
		ret = optimize.fmin(
				fit_functions.fit_GH,
				x0=x0,
				args=(self,sim),
				xtol=self.xtol,
				ftol=self.ftol,
				disp=True)

		print "Fitting of dG and dH parameters concluded:"
		print " dG0 parameters:"
		for (i,f) in enumerate(ret[:len(sim.dG)]):
			print "  %i : %.5E" % (i,f)
		print " dH0 parameters:"
		for (i,f) in enumerate(ret[len(sim.dG):]):
			print "  %i : %.5E" % (i,f)

		if writeback:
			for T in sim.experiments.keys():
				for E in sim.experiments[T]:
					E.dG = ret[:len(sim.dG)]
					E.dH = ret[len(sim.dG):]

		return (ret[:len(sim.dG)], ret[len(sim.dG):])

	def fit_GHC(self, sim, writeback=False):
		"""
		Find optimal dG, dH, and dCp values for all experimental datasets
		"""

		if len(sim.experiments) < 2:
			print "Error:   Cannot fit dCp with less than one experimental temperature"
			return

		x0 = list(sim.dG) + list(sim.dH) + list(sim.dCp)
		ret = optimize.fmin(
			fit_functions.fit_GHC,
			x0=x0,
			args=(self,sim),
			xtol=self.xtol,
			ftol=self.ftol,
			disp=True)

		print "Fitting of dG, dH, and dCp parameters concluded:"
		print " dG0 parameters (@ %.3fK):" % (sim.T_ref)
		for (i,f) in enumerate(ret[:len(sim.dG)]):
			print "  %i : %.5E" % (i,f)
		print " dH0 parameter: (@ %.3fK):" % (sim.T_ref)
		for (i,f) in enumerate(ret[len(sim.dG):2*len(sim.dG)]):
			print "  %i : %.5E" % (i,f)
		print " dCp parameters:"
		for (i,f) in enumerate(ret[2*len(sim.dG):]):
			print "  %i : %.5E" % (i,f)

		if writeback:
			for T in sim.experiments.keys():
				for E in sim.experiments[T]:
					E.dG = ret[:len(sim.dG)]
					E.dH = ret[len(sim.dG):]

		return (ret[:len(sim.dG)], ret[len(sim.dG):2*len(sim.dG)], ret[2*len(sim.dG):])

	def fit_G_opt_H(self, sim, writeback=False):
		(dG,dH) = opt_functions.fit_G_opt_H(self,sim)

		print "Fitting of dG parameters with dH values optimized per-experiment concluded:"
		print " dG0 parameters:"
		for (i,f) in enumerate(dG):
			print "  %i : %.5E" % (i,f)

		print "dH table:"
		for title in dH:
			print "%s :\t" % title,
			for f in dH[title]:
				print "\t%.5E" % f,
			print ""

		if writeback:
			for T in sim.experiments.keys():
				for E in sim.experiments[T]:
					E.dG = dG
					E.dH = dH[ E.title ]

		return(dG,dH)








