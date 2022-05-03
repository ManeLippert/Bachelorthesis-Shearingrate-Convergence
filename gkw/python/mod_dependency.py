#!/usr/bin/env python

# The function below is based on code found on
# http://fusionplasma.co.uk/2013/10/30/Graphing-Fortran-Projects/
# at 2014/10/07.
import re

def build_graph(files,exclude=['hdf5','h5lt','mpi']):
	"""Build a dot graph of the dependency of the modules in the list of files,
	excluding modules named in the list exclude"""

	# Start the graph
	graph = "digraph G {\n"
	deps = {}

	# Pattern for getting the module name.
	p_mod = re.compile("^(?:module|program) ([a-zA-Z0-9_]*)",
	re.IGNORECASE+re.MULTILINE)
	# Pattern for getting the use statements.
	p_use = re.compile("USE ([a-zA-Z0-9_]*),",
	re.IGNORECASE+re.MULTILINE)

	for mod in files:
		with open(mod,'r') as f:
			contents = f.read()

		# Get the true module name
		# Careful! It only gets the first module in the file!
		#if hasattr(re.search(p, contents), 'group'):
		mod_name = re.search(p_mod, contents).group(1)
		#else:
			#print contents

		# Find all the USE statements and get the module
		deps[mod_name] = re.findall(p_use, contents)

	for k in deps:
		# Remove duplicates and unwanted modules
		deps[k] = list(set(deps[k]) - set(exclude))
		for v in deps[k]:
			# Add the edges to the graph
			edge = k + " -> " + v + ";\n"
			graph = graph + edge
 
	# Close the graph and return it
	graph = graph + "}"
	return graph

def build_text_array(files,exclude=['hdf5','h5lt','mpi']):
	"""Build a text array of the dependency of the modules in the list of files,
	excluding modules named in the list exclude"""

	deps = {}

	# Pattern for getting the module name.
	p_mod = re.compile("^(?:module|program) ([a-zA-Z0-9_]*)",
	re.IGNORECASE+re.MULTILINE)
	# Pattern for getting the use statements.
	p_use = re.compile("USE ([a-zA-Z0-9_]*),",
	re.IGNORECASE+re.MULTILINE)

	for mod in files:
		with open(mod,'r') as f:
			contents = f.read()

		# Get the true module name
		# Careful! It only gets the first module in the file!
		mod_name = re.search(p_mod, contents).group(1)

		# Find all the USE statements and get the module
		deps[mod_name] = re.findall(p_use, contents)

	index = 1
	text_array = ""
	head_list = []
	for k in deps:
		# Remove duplicates and unwanted modules
		deps[k] = list(set(deps[k]) - set(exclude))
		text_array = text_array + '%(number)4d ' % {"number": index}
		for v in deps:
			if v in deps[k]:
				text_array = text_array + "x"
			else:
				text_array = text_array + " "
		text_array = text_array + "\n"
		head_list.append('%(number)4d' % {"number": index})
		index = index + 1

	text_array = text_array + "\n"
	index = 1

	# Loop for the key.
	for k in deps:
		text_array = text_array + '%(number)4d' % {"number": index} + " = " + k + "\n"
		index = index + 1

	text_array = text_array + "Note: Read column numbers vertically.\n"
	text_array = text_array + "Module in row, depends on modules in columns with 'x'.\n"

	header = ""
	# Double loop for creating the header.
	for i in range(1,4):
		# Add the space for the row numbering.
		header = header + "     "
		for k in head_list:
			header = header + k[i]
		header = header + "\n"
	# Additional line to sepe
	header = header + "\n"

	text_array = header + text_array

	return text_array

# Wanted: Text, at the beginning tabular structure. on top and on left are
# numbers (both with 4 places). The ones on top written vertically, thus each
# position can be described by two numbers. A 'x' at a certain location means
# that the module defined by the number on the left depends on the module
# defined by the number on the top, directly by a use statement. Below the
# tabular structure, there is a listing in the form number = modulename, thus
# allowing one to make the relation.
#
# for the array part: Loop over 'deps' (= rows), inside another loop over
# 'deps' (columns), if present in deps[k] ad 'x' else ' '.
#
# Missing: Numbers for columns.
# building list with the numbers, then iterating with switched dimensions?

if __name__ == '__main__':
	from sys import argv
	graph = build_graph(argv[1:-1], [])
	#graph = build_text_array(argv[1:-1], [])

	# Print the graph to the standard-output.
	print(graph)
