from itcsimlib.thermo import *
from itcsimlib.model_drakon import *

# Autogenerated with DRAKON Editor 1.29
class Model(DRAKONIsingModel):
    description = "This is an Ising-based model with energies defined by a DRAKON graphical flow diagram."

    def setup(self):
        #item 136
        self.initialize(nsites=8, circular=True)
        #item 30
        self.add_parameter("dG_bind", type="dG")
        #item 32
        self.add_parameter("dH_bind", type="dH")
        #item 35
        self.add_parameter("dCp_bind", type="dCp")


    def site(self, config, site):
        #item 97
        if self.occupied(config, site):
            #item 90
            self.add_dG(config, 'dG_bind', dH='dH_bind', dCp='dCp_bind')
            #item 143
            self.add_dH(config, 'dH_bind', dCp='dCp_bind')
        else:
            pass


