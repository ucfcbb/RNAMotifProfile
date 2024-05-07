def read_isosteric_scoring_data(input_fname):
	rows, cols = (64, 64)
	isosteric_substitution_matrix = [] #[[0 for i in range(cols)] for j in range(rows)]

	fp = open(input_fname)
	lines = fp.readlines()
	fp.close()

	for line in lines:
		if line.startswith('#'):
			continue
		# pieces = line.strip().split('\t')
		# pieces = list(map(float, line.strip().split('\t')))
		pieces = list(map(lambda x: float(x) * 100.0, line.strip().split('\t')))
		isosteric_substitution_matrix.append(pieces)

	return isosteric_substitution_matrix

def read_stacking_scoring_data(input_fname):
	rows, cols = (4, 4)
	stacking_substitution_matrix = [] #[[0 for i in range(cols)] for j in range(rows)]

	fp = open(input_fname)
	lines = fp.readlines()
	fp.close()

	for line in lines:
		if line.startswith('#'):
			continue
		# pieces = line.strip().split('\t')
		# pieces = list(map(float, line.strip().split('\t')))
		pieces = list(map(lambda x: float(x) * 100.0, line.strip().split('\t')))
		stacking_substitution_matrix.append(pieces)

	return stacking_substitution_matrix

def read_nucleotide_scoring_data(input_fname):
	rows, cols = (4, 4)
	nucleotide_substitution_matrix = [] #[[0 for i in range(cols)] for j in range(rows)]

	fp = open(input_fname)
	lines = fp.readlines()
	fp.close()

	for line in lines:
		if line.startswith('#'):
			continue
		# pieces = line.strip().split('\t')
		# pieces = list(map(float, line.strip().split('\t')))
		# pieces = list(map(lambda x: int(x * 100), pieces))
		pieces = list(map(lambda x: float(x) * 100.0, line.strip().split('\t')))
		nucleotide_substitution_matrix.append(pieces)

	return nucleotide_substitution_matrix