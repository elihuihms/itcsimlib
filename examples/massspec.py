#!/usr/bin/env python

#
# This script demonstrates how an itcsimlib model can be adapted to fit mass spectrometric data
#

from itcsimlib import *
from itcsimlib.model_ising import NonAdditive
from itcsimlib.mass_spec import MSExperiment, MSModel

itc_model = NonAdditive(nsites=11, circular=1, units="J")

ms_model = MSModel(itc_model)
print ms_model

ms_model.set_params(dGX=-27000,dGY=-27000,dGZ=-30000)

sim = ITCSim(verbose=True, threads=1, units="J")

sim.set_model(ms_model)

exp = MSExperiment('data/TRAP_populations_EDDA.txt')
print exp

sim.add_experiment( exp )

sim.run()

sim.make_plots(hardcopy=True,hardcopytype='png')

sim.experiments[0].export_to_file("massspec.fit")

sim.done()
