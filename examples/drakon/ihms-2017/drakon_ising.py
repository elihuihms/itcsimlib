from itcsimlib.thermo import *
from itcsimlib.model_drakon import *

# Autogenerated with DRAKON Editor 1.29
class Model(DRAKONIsingModel):
    description = "This is an Ising-based model with energies defined by a DRAKON graphical flow diagram."

    def setup(self):
        #item 111
        num_sites = 11
        circular = True
        #item 136
        self.initialize(nsites=num_sites,circular=circular)
        #item 30
        self.add_parameter("dG_0",type="dG")
        self.add_parameter("dG_oe",type="dG")
        self.add_parameter("dG_oo",type="dG")
        #item 32
        self.add_parameter("dH_0",type="dH")
        self.add_parameter("dH_oe",type="dH")
        self.add_parameter("dH_oo",type="dH")
        #item 35
        self.add_parameter("dCp_0",type="dCp")
        self.add_parameter("dCp_oe",type="dCp")
        self.add_parameter("dCp_oo",type="dCp")


    def site(self, i, j):
        #item 97
        if self.occupied(i,j) == True:
            #item 90
            self.add_dG(i, 'dG_0', dH='dH_0', dCp='dCp_0')
            self.add_dH(i, 'dH_0', dCp='dCp_0')
            #item 95
            if self.occupied(i,j-1) == True:
                #item 91
                if self.occupied(i,j+1) == True:
                    #item 108
                    self.add_dG(i, 'dG_oo', dH='dH_oo', dCp='dCp_oo')
                    self.add_dH(i, 'dH_oo', dCp='dCp_oo')
                else:
                    #item 109
                    self.add_dG(i, 'dG_oe', dH='dH_oe', dCp='dCp_oe')
                    self.add_dH(i, 'dH_oe', dCp='dCp_oe')
                    #item 99
                    if self.occupied(i,j+1) == True:
                        #item 110
                        self.add_dG(i, 'dG_oe', dH='dH_oe', dCp='dCp_oe')
                        self.add_dH(i, 'dH_oe', dCp='dCp_oe')
                    else:
                        pass
            else:
                #item 99
                if self.occupied(i,j+1) == True:
                    #item 110
                    self.add_dG(i, 'dG_oe', dH='dH_oe', dCp='dCp_oe')
                    self.add_dH(i, 'dH_oe', dCp='dCp_oe')
                else:
                    pass
        else:
            pass


