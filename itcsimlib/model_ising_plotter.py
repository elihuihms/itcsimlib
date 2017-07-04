from math import pi,sin,cos
from pyx import *

def draw_ising_lattices(path, model, size=1.0):
	c = canvas.canvas()
	model.set_energies(273.15,273.15) # ensure that this is run at least once to populate config_terms
	
	for i in xrange(model.nsites+1):
		configurations = {} # keyed by energy, (key,weight)
		for j in xrange(model.nconfigs):
			if model.bound[j] == i:
				if model.gibbs[j] in configurations.keys():
					configurations[ model.gibbs[j] ][1]+=1
				else:
					configurations[ model.gibbs[j] ] = [model.configs[j],1]
		
		for j,energy in enumerate(sorted(configurations.keys())):
			if model.circular:
				_draw_circular_lattice(c,
					(2.2*size*j,0.0-2.2*size*i),size,model.nsites,
					key=configurations[energy][0],energy="%.1f"%(energy),degeneracy=configurations[energy][1])
			else:
				_draw_linear_lattice(c,
					(1.1*size*j,1.1*size*i*model.nsites),size,model.nsites,
					key=configurations[energy][0],energy="%.1f"%(energy),degeneracy=configurations[energy][1])
					
	c.writePDFfile(path)

def _draw_linear_lattice(c,loc,w,n,key=None,energy=None,degeneracy=None):	
	x,y = loc
	sx,sy = w,n*w

	x0,x1 = x-(sx/2.0),x+(sx/2.0)
	y0,y1 = y-(sy/2.0),y+(sy/2.0)
	rect = box.polygon([(x0, y0),(x1, y0),(x1, y1),(x0, y1)])
	c.stroke(rect.path(), [style.linewidth.Thick]) # deformer.smoothed(radius=w*2)])

	for i in xrange(0,n):
		if key[i] > 0:
			c.fill(path.circle(x,y+(i*w)-(sy/2.0)+(w/2.0),w*0.3), [color.rgb.blue])
		c.stroke(path.circle(x,y+(i*w)-(sy/2.0)+(w/2.0),w*0.4), [style.linewidth.Thick])

def _draw_circular_lattice(c,loc,r,n,key=None,energy=None,degeneracy=None):
	x,y = loc
	chord = 2.0*r*sin(pi/n)
	site_r = (r-(chord/2.0))*sin(pi/n)
	outer_r = site_r / sin(pi/n) + site_r
	inner_r = site_r / sin(pi/n) - site_r

	c.stroke(path.circle(x,y,outer_r), [style.linewidth.Thick])
	c.stroke(path.circle(x,y,inner_r), [style.linewidth.Thick])

	if energy != None:
		texrun = text.defaulttexrunner
		t = texrun.text(x, y, r"%s"%str(energy), [text.halign.boxcenter, text.valign.middle] )
		if unit.tocm(t.width) < 2.0*inner_r:
			c.insert(t)
	
	if degeneracy != None:
		t = c.text(x+(r*0.68),y-(r*0.68), r"x%s"%str(degeneracy), [text.halign.boxleft, text.valign.top])

	for i in xrange(0,n):
		circ_x = x+(cos(2.0*pi/n*i)*(r-(0.5*chord)))
		circ_y = y+(sin(2.0*pi/n*i)*(r-(0.5*chord)))
		if key[i] > 0:
			c.fill(path.circle(circ_x,circ_y,site_r*0.6), [color.rgb.blue])
		c.stroke(path.circle(circ_x,circ_y,site_r*0.8), [style.linewidth.Thick])
