#!/bin/bash

# --------------------------------
# 1. commandline arguments
# --------------------------------
timelimit=1				    	# Time limit of simulation in hours
seed=1234						# random seed to generate initial conditions 
beg_idx=0 						# the index at which the configurations start to be loaded.
Verbosity=3						# verbosity of output messages  

# Interactive arguments.
fname=test		# save file

# --------------------------------
# 2. Parameters of pair potential.
# --------------------------------
d=1                            # space dimension
K1=0.0                         # lower bound on the exclusion region (in unit number density, set it 0.0 for stealthy hyperuniform)
chi=0.3				# stealthiness parameter (assume stealthy hyperuniformity)
S0=0.				# S0 > 0: equiluminous. S0 = 0.0: Stealthy
vareps0=100.0			# relative strength of the soft-core repulsion. (0 means no soft-core repulsions.)
phi_fake=0.15                      # (Do not touch) fictitious packing fraction; the particle radius is computed from this value in the unit number density
sigma=0.20                     # Exclusion radius of the soft-core repulsion in unit number density; this value must be larger than the particle diameter

# --------------------------------
# 3. Computational parameters
# --------------------------------
threads=2                      # Number of threads in openmp
N=100	                       # Number of particles
Nc=10                          # Number of configurations
initial=random                 # Choice of initial condition (random/input)

# Minimization parameters
# - These parameters can be skipped. They are also insensitive in order and case-insensitive.
eps0="tolerance 1e-14"		# energy tolerance for ground states in the unit of v0. Default value is 1e-14.
max_eval="maxsteps 100000" 	# the maximum number of evaluations before quitting calculations. Default value is 10000.
algorithm="algorithm LBFGS"	# Minimization algorithm. Default option is LBFGS. Other available options: LocalGradientDescent, ConjugateGradient, SteepestDescent, MINOP


# ------------- do not touch ---------
#K2=0.2 	                   # upper bound on the exclusion region (in unit number density)
pi=`echo "4*a(1)" | bc -l`
if [ "$d" == 1 ]; then
v=2.0
elif [ "$d" == 2 ]; then
v=${pi} 
elif [ "$d" == 3 ]; then 
v=`echo "4.*${pi}/3."| bc -l `
fi 
K2=`echo "2*${pi}*e(l(2*${d}*${chi}/${v})/${d})" | bc -l `
a=`echo "e(l(${phi_fake}/${v})/${d})" | bc -l`
K1a=`echo "${K1}*$a"| bc -l `
K2a=`echo "${K2}*$a"| bc -l `
# ------------------------------------


# Filenames
#exe=${HOME}/code/CCO/CCOptimization/soft_core_stealthy2.out
exe=./CCO.out

echo "# Run the following command:"
#echo "${exe} ${timelimit} ${eps0} ${max_eval} <<< \"${d} 0 ${Ka} 0.0 1.0 ${sigma} ${phi2} ${threads} ${N} ${Nc} ${Nc_MD} ${initial} ${fname} ground\""

# ground states from Random Initial Conditions
${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${K1a} ${K2a} $S0 $vareps0 ${sigma} ${phi_fake} \
${threads} ${N} ${Nc} random ${fname}_GS \
ground $eps0 $max_eval $algorithm run" 


# ground states from Input Initial Conditions
#loadconfig=${fname}_GS # the name of loaded ConfigPack file.
#K1a=0.				   # change in K1 value.
#${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${K1a} ${K2a} $S0 $vareps0 ${sigma} ${phi_fake} \
#${threads} ${N} ${Nc} input ${loadconfig} ${fname}_GS2 \
#ground $eps0 $max_eval $algorithm run" > log3

