# Summary
	C++ implementations of computing statistics for point patterns.
	Most of the codes were written by Ge Zhang, a former graduate student in the Torquato group.

# Requirements

	* gsl library
	* boost library
		* The `boost/random/*` must be available
		
	* `../core` from git@github.com:jaeukk/cores.git 
	* Make

# Compilation
	>> export core=/path/to/core/from/git

	Then, the code can be compiled with the provided makefile using the standard `make` command.
	This creates an executable ~/EXC/PairStat.out , which is a CLI for computing some predefined statistics.

# Sample Execution & Output

	If run without command line arguments, using

	

	