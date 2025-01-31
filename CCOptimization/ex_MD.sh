#!/bin/bash

# commandline arguments
timelimit=1				    	# Time limit of simulation in hours
seed=1234						# random seed to generate initial conditions 
beg_idx=0 						# the index at which the configurations start to be loaded.
Verbosity=3						# verbosity of output messages  

# Interactive arguments.
# Parameters of pair potential.
d=1                            # space dimension
K1=0.01                    # lower bound on the exclusion region (in unit number density)
K2=0.047124 	                   # upper bound on the exclusion region (in unit number density)
S0=0.						   # S0 > 0: equiluminous. S0 = 0.0: Stealthy
vareps0=1.0							# relative energy of the soft-core repulsion. (0 means no soft-core repulsions.)
phi_fake=0.15                      # Packing fraction; the particle radius is computed from this value in the unit number density
sigma=0.20                     # Exclusion radius in unit number density; this value must be larger than the particle diameter

# Computational parameters
threads=2                      # Number of threads in openmp
N=100	                         # Number of particles
Nc=10                          # Number of configurations
initial=random                 # Choice of initial condition (random/input)

# MD simulation parameters 
# - These parameters can be skipped. They are also insensitive in order and case-insensitive.
T="temperature 1e-3" 		# set dimensionless temperature for MD simulations. -1 is the default value, implying the sufficiently small temperature. Otherwise, it must be positive.
dt="timestep 0.05"			# set dimensionless time step. 0.01 is the default value.
int_samp="samplesteps 10000" # set sample interval along the MD trajectory. 1e5 is the default value. 
num_equi="numequilibration 100" # set the number of steps in MD simulations spent for equilibration. Those configurations are discarded. If you want to turn on "AutoTimeStep", this quantity must be positive!!
auto_dt="autotimestep"		# make the code automatically determine the MD time step. By defualt, it uses a fixed time step.  

# ------------------ Do not touch --------------
if [ "$d" == 1 ]; then
v=2.0
elif [ "$d" == 2 ]; then
v=3.1415926535898
elif [ "$d" == 3 ]; then 
v=`echo "4.*${v}/3."| bc -l `
fi 
a=`echo "e(l(${phi_fake}/${v})/${d})" | bc -l`
K1a=`echo "${K1}*$a"| bc -l `
K2a=`echo "${K2}*$a"| bc -l `
# -----------------------------------------------

# Filenames
#exe=${HOME}/code/CCO/CCOptimization/soft_core_stealthy2.out
exe=./CCO.out
fname=test

echo "# Run the following command:"
#echo "${exe} ${timelimit} ${eps0} ${max_eval} <<< \"${d} 0 ${Ka} 0.0 1.0 ${sigma} ${phi2} ${threads} ${N} ${Nc} ${Nc_MD} ${initial} ${fname} ground\""

# MD simulations at temperature $T from random initial conditions.
${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity}  <<< "${d} ${K1a} ${K2a} $S0 $vareps0 ${sigma} ${phi_fake} \
${threads} ${N} ${Nc} random ${fname}_MD \
MD $T $dt $int_samp $num_equi $auto_dt run" > log_MD

