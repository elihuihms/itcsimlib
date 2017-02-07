import math

from Tkinter	import *

class Manipulator(Tk):

	def __init__(self, sim, experiments=None, params=None, callback=None, title="Value Manipulator"):
		Tk.__init__(self)
		self.protocol('WM_DELETE_WINDOW', self.close)

		self.visualizer = Visualizer(self, title="Plot Window")

		self.sim = sim
		self.callback = callback
		self.title(title)
		
		if experiments:
			self.experiments = experiments
		else:
			self.experiments = self.sim.get_experiments()
		if params:
			self.params = params
		else:
			self.params = self.sim.get_model_params()

		for e in self.experiments:
			self.visualizer.add_experiment(e)
			
		self.sim.run(self.experiments)
		self.visualizer.make_figures()

		self.size = 0
		self.sliders = []
		self.entries = []

		grid_counter = 0
		for p in self.params:
			value = self.sim.get_model_param(p)
			self.sliders.append( Scale(self, from_=-1.0*math.fabs(value), to=2.0*math.fabs(value), label=p, length=500, orient=HORIZONTAL) )
			self.sliders[-1].grid( row=self.size, column=1 )
			self.sliders[-1].set( value )

			self.entries.append( Entry(self) )
			self.entries[-1].grid( row=self.size, column=3 )
			self.entries[-1].insert(0,self.sliders[-1].get())

			self.size+=1

		self.bind('<B1-Motion>',self.slide_update)
		self.bind('<ButtonRelease-1>',self.slide_rescale)
		self.bind('<Return>',self.entry_update)

		for i in xrange(self.size):
			value = float(self.entries[i].get())
			self.sliders[i].config(to	=1.5*value)
			self.sliders[i].config(from_=0.5*value)
			self.sliders[i].set(value)

		self.mainloop()

	def _update_figures(self):
		if self.callback != None:
			self.callback([float(self.entries[i].get()) for i in xrange(self.size)])
		else:
			for i,p in enumerate(self.params):
				self.sim.set_model_param(p,float(self.entries[i].get()))

		self.sim.run(self.experiments)
		self.visualizer.update_figures()

	def slide_update(self, event):

		if event.widget in self.sliders:
			index = self.sliders.index(event.widget)
			value = event.widget.get()
			self.entries[index].delete(0,END)
			self.entries[index].insert(0,value)

			event.widget.set(value)
			self._update_figures()

	def slide_rescale(self,event):

		if event.widget in self.sliders:
			value = event.widget.get()
			scale = event.widget.cget("to") -event.widget.cget("from")

			if value == event.widget.cget("to"):
				event.widget.config(to		=value+scale)
				event.widget.config(from_	=value)
			elif value == event.widget.cget("from"):
				event.widget.config(to		=value)
				event.widget.config(from_	=value-scale)

	def entry_update(self, event):

		index = self.entries.index(event.widget)
		value = float(event.widget.get())

		self.sliders[index].config(to	=1.5*value)
		self.sliders[index].config(from_=0.5*value)
		self.sliders[index].set(value)

		self._update_figures()
		
	def close(self):
		self.destroy()
		self.quit()

class Visualizer(Toplevel):

	def __init__(self,master,title):
		Toplevel.__init__(self,master)
		self.master = master

		from __init__ import MATPLOTLIB_BACKEND
		if MATPLOTLIB_BACKEND != None:
			print "manipulator: Setting matplotlib backend to \"TkAgg\"."
			
		import matplotlib
		matplotlib.use("TkAgg")

		from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
		from matplotlib.figure import Figure
		from matplotlib import pyplot

		self.title(title)
		self.resizable(True,True)
		self.fig = pyplot.figure()
		pyplot.ion()

		self.canvas = FigureCanvasTkAgg(self.fig, master=self)
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
		self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
		self.update()

		self.experiments = []

	def add_experiment(self,experiment):
		self.experiments.append(experiment)

	def make_figures(self):
		size = math.ceil(math.sqrt(len(self.experiments)))
		self.axs,self.exp,self.fit = [],[],[]
		for i,e in enumerate(self.experiments):
			self.axs.append( self.fig.add_subplot(size,size,i+1) )
			ratios = [ e.Concentrations[i][e.syringeRef]/e.Concentrations[i][e.cellRef] for i in xrange(e.npoints) if i not in e.skip ]
			self.exp.append( self.axs[-1].scatter(ratios,e.dQ_exp,c='#000000') )
			self.fit.append( self.axs[-1].plot(ratios,e.dQ_fit,c='#FF0000')[0] )
		self.canvas.draw()

	def update_figures(self):
		for i,e in enumerate(self.experiments):
			self.fit[i].set_ydata(e.dQ_fit)
		self.canvas.draw()


