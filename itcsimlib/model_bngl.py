import os
import tempfile
import shutil
import pexpect
import xml.etree.ElementTree as ET

from collections import OrderedDict

from .thermo	import *
from .itc_model	import ITCModel

	
class BNGL(ITCModel):

	def __init__(self,model,interpreter):
		ITCModel.__init__(self)

		self.model = model
		self.interpreter = interpreter

		# need to briefly fire up the BNG interpreter in order to generate the XML file
		self.start()
		self.parse_network_xml()
		self.stop()

		for p in self.bngl_parameters:
			if p[:2] == 'E_': # crude check for energy terms, is there a better way to do this?
				self.add_parameter( "dG_%s"%(p[2:]),	'dG',	default=self.bngl_parameters[p],description='' )
				self.add_parameter( "dH_%s"%(p[2:]),	'dH',	description='' )
				self.add_parameter( "dCp_%s"%(p[2:]),	'dCp',	description='' )
			elif p[:2] == 'k_': # crude check for rate terms
				self.add_parameter( p,	'k',	default=self.bngl_parameters[p],description='' )
				
	def start(self):
		# set up the BNG environment
		self.workdir = tempfile.mkdtemp()
		self.basepath = os.path.join(self.workdir,os.path.splitext(os.path.basename(self.model))[0])
		
		shutil.copyfile(model,"%s.bngl"%self.basepath)
		self.log = open("%s.log"%self.basepath,'w+')		

		# start the BNGL interpreter
		self.start_bng()		
					
	def stop(self):
		self.bng.sendline("action quit()")
		self.bng.close()
		self.bng = None

		self.log.close()
			
		shutil.rmtree( self.workdir )
		
	def dump_log(self,file=None):
		self.log.seek(0)
		if file:
			f = open(file)
			f.write(self.log.read())
			f.close()
		else:
			print(self.log.read())

	def send_bng(self,cmd):
		self.bng.sendline(cmd)
		self.bng.expect("BNG>")

	def start_bng(self):
		self.bng = pexpect.spawn("%s --console"%self.interpreter,cwd=self.workdir,logfile=self.log)
		self.bng.expect("BNG>")
		self.send_bng("load \"%s.bngl\""%self.basepath)
		self.send_bng("action generate_network({overwrite=>1})")
		self.send_bng("action writeXML()")
		
	def parse_network_xml(self):
		tree	= ET.parse( "%s.xml"%self.basepath )
		root	= tree.getroot()
		model	= root[0]

		namespace = root.tag[:-4]
		def mktag(s):
			return "%s%s"%(namespace,s)

		self.bngl_parameters = {}
		for p in model.find(mktag("ListOfParameters")):
			if p.attrib['type'] == 'Constant':
				self.bngl_parameters[ p.attrib['id'] ] = float(p.attrib['value'])
			else:
				self.bngl_parameters[ p.attrib['id'] ] = p.attrib['value']
			
		for m in model.find(mktag("ListOfMoleculeTypes")):
			tmp = [c.attrib['id'] for c in m.find(mktag("ListOfComponentTypes"))]
			self.add_component("%s(%s)"%(m.attrib['id'],','.join(tmp)))

		self.model_species = OrderedDict()
		for s in model.find(mktag("ListOfSpecies")):
			self.model_species[ s.attrib['name'] ] = {}
			self.model_species[ s.attrib['name'] ]['id'] = s.attrib['id']
	
			energies = s.find(mktag("ListOfEnergies"))
			if energies:
				self.model_species[ s.attrib['name'] ]['energy_expressions'] = [e.attrib['expression'] for e in energies]
			else:
				self.model_species[ s.attrib['name'] ]['energy_expressions'] = []
				
	def parse_cdat(self):
		try:
			h = open( "%s.cdat"%self.basepath )
		except e as Exception:
			self.dump_log()
			raise e

		ret = []				
		for l in h.readlines()[1:]:
			ret.append( map(float,l.split()[2:]) ) # column 0 = empty, column 1 = time
		h.close()
		return ret
		
	def Q(self,T_ref,T,concentrations,steptime=1E3,steppoints=1):			
		self.send_bng("action resetConcentrations()")
		
		dH = {}
		for p in self.params:
			if self.get_param_type(p) == 'k': # rate
				self.send_bng("action setParameter(\"%s\",%.5E)"%(p,self.params[p]))				
			if self.get_param_type(p) == 'dG': # free energy
				self.send_bng("action setParameter(\"E_%s\",%.5E)"%(p[3:],dG[p]))
			elif self.get_param_type(p) == 'dH': # enthalpy
				dH[p] = dH_vant_Hoff( self.params[p], self.params["dCp_%s"%p[3:]], T, T_ref )
				
		enthalpies,prev_concentrations = {},{}
		for name in self.model_species:
			enthalpies[name] = sum( eval(e,self.bngl_parameters) for e in self.model_species[s]['energy_expressions'] )				
			prev_concentrations[s] = 0.0,0.0 # total conc,free conc
		
		Q,counter = [0.0]*(len(concentrations)*steppoints),0
		for i,concs in enumerate(concentrations):
			
			for name,value in concs.items():
				self.send_bng("action setConcentration(\"%s\",%.5E)"%(name,(value - prev_concentrations[name][0]) +prev_concentrations[name][1]))
				
				#increase total conc
				prev_concentrations[name][0]+=value
			
			self.send_bng("action simulate({method=>\"ode\",t_start=0,t_end=>%i,t_steps=>%i})"%(steptime,steppoints))
			timecourse = self.parse_cdat()
			
			for j,name in enumerate(self.model_species):
			
				# increment the heat content in the cell from the enthalpy of each component
				for k in range(steppoints):
					Q[counter+k] += timecourse[k][j] * enthalpies[name]
				
				#set final free conc
				prev_concentrations[name][1] = timecourse[k][j]
			
			counter+=steppoints
		return Q
				