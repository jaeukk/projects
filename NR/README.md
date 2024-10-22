# Summary
A code to compute a probability distribution $P[N(R)]$ associated with the number of points $N(R)$ within a spherical window of radius $R$ for an ensemble of point patterns.
It stores all instances of $N(R)$ for a given value of $R$ in ".bin" files.


Please cite the following references:
- S. Torquato, J. Kim, and M. A. Klatt, Local Number Fluctuations in Hyperuniform and Nonhyperuniform Systems: Higher-Order Moments and Distribution Functions, Physical Review X, 11, 021028 (2021).
- M. A. Klatt, J. Kim, T. E. Gartner III, and S. Torquato, Local number fluctuations in ordered and disordered phases of water across temperatures: Higher-order moments and degrees of tetrahedrality , The Journal of Chemical Physics, 160 204502 (2024).


# Requirements

	* gsl library
	* boost library
		* The `boost/random/*` must be available
		
	* `../core` from git@github.com:jaeukk/cores.git 
	* Make

# Compilation
	>> export core=/path/to/core/from/git

	Then, the code can be compiled with the provided makefile using the standard `make` command.
	This creates an executable ./NR.out , which is a CLI for computing number distributions.

# Sample Execution & Output

	PrintPDF=1			# 1=print probability distribution functions; do not print otherwise.
	PrintCDF=1			# 1=print cumulative probability distribution functions; do not print otherwise.
	ImplementPBC=1		# 1=implement periodic boundary conditions; do not implement PBC otherwise.
	SamplingOption=-1	# <0 = use uniform sampling windows; Use windows on a regular cubic grid the specified spacing, otherwise.

	Name=<Name>			# Name of a ConfigPack file name to load
	Idx_beg=0			# The index of the first configuration.
	Nc=10				# The number of configurations after Idx_beg
	seed=12312			# Random seed
	N_window=100		# The number of windows
	Rmax=30				# The maximum radius of windows.
	dR=0.1				# Increment in window radius
	threads=10			# Number of threads in OpenMP
	binfilename=<name2>	# The name of temporary binary (.bin) files
	output=<name3>		# The name of output files (e.g., distribution functions, moments, L2 norms, ....).

	./NR.out ${PrintPDF} ${PrintCDF} ${ImplementPBC} ${SamplingOption} <<< \ 
	"$Name 
	$Idx_beg	
	$Nc			
	$seed		
	$N_window	
	$Rmax		
	$dR			
	$threads	
	$binfilename
	$output	"
