# Summary
	C++ implementations of managing the ConfigPack files to store a list of point configurations under the periodic boundary conditions.
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
	This creates an executable ~/EXC/manager.out , which is a CLI for handling ConfigPack files.

# Sample Execution & Output

	>> ./manager.out $option [other arguments]

	* Check very simple statistics in the ConfigPack file
	option=0 # 
	name=test #name of the ConfigPack file (without the file extension)
	ex) ./manager.out $option $name

	* Convert a portion of configuration in the ConfigPack file into the txt format
	option=1 
	name=test # name of the ConfigPack file (without the file extension)
	idx_beg=0 #index of the starting configuration
	idx_end=10 #the number of configurations
	ex) ./manager.out $option $name $idx_beg #idx_end

	* Add a list of txt files to a ConfigPack file.
	option=2
	name1=test_%02d  #the name of txt files without the extension
	idx_beg=0 # the starting index
	Nc=10 # the number of configurations
	savename=output #the name of the ConfigPack file (without the extension ConfigPack)
	ex) ./mananger.out $option $name1 $idx_beg $Nc $savename

	In option 2, if you want to convert a list of text files whose name is file_01.txt, ... file_09.txt, the usage will be
	>>./mananger.out 2 ./file_%02d 1 9 [the name of the ConfigPack file]


	

	