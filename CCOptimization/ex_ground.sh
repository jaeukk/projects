#!/bin/bash

# commandline arguments
timelimit=1				    	# Time limit of simulation in hours
eps0=1e-18                     # upper bound on energy
max_eval=10000000  				# max number of energy evaluations.
seed=1234						# random seed to generate initial conditions
beg_idx=0 						# the index at which the configurations start to be loaded.
Verbosity=3						# verbosity of output messages
algorithm=BFGS					# numerical minimization algorithm. BFGS is by default
T=1e-3 							# dimensionless temperature for MD simulations. -1 is the default value, implying the sufficiently small temperature. Otherwise, it must be positive.


# Interactive arguments.
d=1                            # space dimension
K1=0.047124                    # lower bound on the exclusion region (in unit number density)
K2=0.01 	                   # upper bound on the exclusion region (in unit number density)
S0=0.						   # S0 > 0: equiluminous. S0 = 0.0: Stealthy
v0=1.0							# relative energy of the soft-core repulsion. (0 means no soft-core repulsions.)
phi_fake=0.15                      # Packing fraction; the particle radius is computed from this value in the unit number density
sigma=0.20                     # Exclusion radius in unit number density; this value must be larger than the particle diameter

threads=8                      # Number of threads in openmp
N=100	                         # Number of particles
Nc=10                          # Number of configurations
Nc_MD=0                        # number of steps in MD simulations
initial=random                 # Choice of initial condition (random/input)

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

# Filenames
#exe=${HOME}/code/CCO/CCOptimization/soft_core_stealthy2.out
exe=./CCO.out
fname=test

echo "# Run the following command:"
#echo "${exe} ${timelimit} ${eps0} ${max_eval} <<< \"${d} 0 ${Ka} 0.0 1.0 ${sigma} ${phi2} ${threads} ${N} ${Nc} ${Nc_MD} ${initial} ${fname} ground\""

# MD simulations at temperature $T from random initial conditions.
${exe} ${timelimit} ${eps0} ${max_eval} ${seed} ${beg_idx} ${Verbosity} ${algorithm} $T <<< "${d} ${K1a} ${K2a} $S0 $v0 ${sigma} ${phi_fake} ${threads} ${N} ${Nc} ${Nc_MD} random ${fname}1 MD" > log1

# ground states from Input Initial Conditions
${exe} ${timelimit} ${eps0} ${max_eval} ${seed} ${beg_idx} ${Verbosity} ${algorithm} <<< "${d} ${K1a} ${K2a} $S0 $v0 ${sigma} ${phi_fake} ${threads} ${N} ${Nc} ${Nc_MD} input ${fname}1 ${fname}2 ground" > log2

# ground states from Random Initial Conditions
${exe} ${timelimit} ${eps0} ${max_eval} ${seed} ${beg_idx} ${Verbosity} ${algorithm} <<< "${d} ${K1a} ${K2a} $S0 $v0 ${sigma} ${phi_fake} ${threads} ${N} ${Nc} ${Nc_MD} random ${fname}3 ground" > log3

