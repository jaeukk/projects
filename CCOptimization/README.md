# Summary
    A C++ code for performing simulations with stealthy pair potential + soft-core repulsion.
    Specifically, it can 1) find numerical ground states, or 2) perform molecular dynamic simulations.

# Requirements
    * gsl library
    * nlopt library (https://nlopt.readthedocs.io/en/latest/NLopt_Installation/)
    * `../cores/`       (inside jaeukk/projects/)  
    * `../PairStat/`    (inside jaeukk/projects/)  
    * `../Potential/`   (inside jaeukk/projects/)  
    * Make

# Compilation
In the makefile, set the paths of three directories (lines 45-47): 
core = /path/to/cores
pair = /path/to/PairStat
pot = /path/to/Potential

Also re-define the variables (IDIR, ODIR, ..., GSLLIB, NLOPTLIB) in lines 49-60.

Then, the code can be compiled with the provided makefile using the standard `make` command.
This creates an executable ${EXEdir}/CCO.out, which is a CLI.

# Usages
1. Finding numerical ground states (i.e., stealthy hyperuniform configurations); see ./ex_ground.sh for example.

    * References (without soft-core repulsions):
        - U. Uche, F. H. Stillinger, and S. Torquato, Phys. Rev. E 70, 046122 (2004).
        - R. D. Batten, F. H. Stillinger, and S. Torquato, J. Appl. Phys. 104, 033504 (2008).
        - G. Zhang, F. H. Stillinger, and S. Torquato, Phys. Rev. E 92, 022119 (2015). 

    * References (with soft-core repulsions):
        - S. Torquato and J. Kim, Soft Matter 21(24), 4898–4907 (2025);
        - J. Kim and S. Torquato, J. Chem. Phys. 163, 024902 (2025).


2. Performing molecular dynamic (MD) simulations; see ./ex_MD.sh for example.
    * References (with soft-core repulsions):
        - G. Zhang, F. H. Stillinger, and S. Torquato, Phys. Rev. E 92, 022119 (2015). 
        - S. Torquato, G. Zhang, and F. H. Stillinger, Phys. Rev. X 5, 021020 (2015). 



