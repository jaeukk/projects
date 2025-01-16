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

1. Compute structure factor.
	omp_threads=1
	config="GetConfigsFunction ReadConfigPack ../configs/integer-lattice"
	OutputName="Z1-test"

./PairStat.out <<< "NumThread $omp_threads
	$config
	SkComputation 6.3 100 0.01
	BraggPeaks 
	CircularKMin 6.25
	OutputPrefix $OutputName
	Calculation
	Exit
	"


	./PairStat.out <<< "NumThread 1
	GetConfigsFunction ReadConfigPack ../configs/integer-lattice
	SkComputation 6.3 100 0.001 
	CircularKMin 6.2 
	BraggPeaks
	OutputPrefix ./Z1-test
	Calculation
	Done
	"

	