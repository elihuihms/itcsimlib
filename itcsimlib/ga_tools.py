import scipy

def get_stats( population, full_output=False ):
	n = len(population)
	ret = { \
		'n': n, \
		'min': min(population.member_chisquare[0:n]), \
		'max': max(population.member_chisquare[0:n]), \
		'mean': scipy.mean(population.member_chisquare[0:n]), \
		'median': scipy.median(population.member_chisquare[0:n]), \
		'stdev': scipy.std(population.member_chisquare[0:n]), \
		}

	if full_output:
		return ret,population.member_chisquare[0:n]
	else:
		return ret

def get_relationships( population ):
	model = population.get_model()
	p_names = [] # name of model parameters
	for type in population.member_param_map[0]:
		for (group,value) in population.member_param_map[0][type]:
			p_names.extend(group)
	
	p_types = set( [model.get_param_type(p) for p in p_names] ) # list of unique model parameter types
	p_table = dict.fromkeys(p_types) # NxN table of parameter correlations, by parameter types
	p_index = dict.fromkeys(p_types) # Lookup table of NxN indices by parameter names, by parameter types
	p_inver = dict.fromkeys(p_types) # Order of parameters in the NxN table (inverted p_index)
	
	for type in p_types:
		p_inver[type] = [p for p in model.get_param_names() if model.get_param_type(p)==type]
		n = len(p_inver[type])
		p_table[type] = scipy.zeros((n,n))
		p_index[type] = dict(zip(p_inver[type],xrange(n)))
		
		for i in xrange(len(population)):
			for (group,value) in population.member_param_map[i][type]:
				n = len(group)
				for i in xrange(n):
					for j in xrange(n):
						p_table[type][ p_index[type][group[i]] ][ p_index[type][group[j]] ] += 1
			
	
		p_table[type] /= len(population)
						
	return p_table,p_inver

def print_relationships( p_table, p_inver, format="%0.3F" ):
	ret = ""
	for type in p_table:
		ret += type
		ret += "\t"
		ret += ("\t".join(p_inver[type]))
		ret += "\n"
		for i in xrange(len(p_inver[type])):
			ret += p_inver[type][i]
			ret += "\t"
			ret += ("\t".join([format%f for f in p_table[type][i]]))
			ret += "\n"
		ret += "\n"
	return ret
