# Summary

A simple code to compare vertices in two Local Isomorphism (LI) classes for 2D quasicrystalline tilings. 


# Requirements

	* gsl library
	* boost library
		* The `boost/random/*` must be available
		
	* `../core` from git@github.com:jaeukk/cores.git 
	* Make

# Compilation
	>> export core=/path/to/core/from/git

	Then, the code can be compiled with the provided makefile using the standard `make` command.
	This creates an executable ${EXEdir}/QC.out, which is a CLI for comparing two LI classes.

# Sample Execution & Output
	
	ref=<Prefix of a ConfigPack file>	# Configurations of the reference LI class
	ref_id=0							# ID in the ConfigPack file.
	comp=<Prefix of a ConfigPack file>	# Configurations of the compared LI class
	comp_id=0							# ID in the ConfigPack file.
	num_threads=1						# The number of threads.
	savename=<output>					# Prefix of the output txt file.

	./QC.out ConfigDiff2 <<< "
	$ref $ref_id
	$comp $comp_id
	$num_threads
	$savename
	"
