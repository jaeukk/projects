#!/bin/bash

# commandline arguments
timelimit=1				    	# Time limit of simulation in hours
seed=1234						# random seed to generate initial conditions 
beg_idx=0 						# the index at which the configurations start to be loaded.
Verbosity=3						# verbosity of output messages  

# Interactive arguments.
# Parameters of pair potential.
d=1                            # space dimension
K1=0.0                    # lower bound on the exclusion region (in unit number density)
chi=0.47 	                   # upper bound on the exclusion region (in unit number density)
S0=0.						   # S0 > 0: equiluminous. S0 = 0.0: Stealthy
vareps0=1.0							# relative energy of the soft-core repulsion. (0 means no soft-core repulsions.)
#phi_fake=0.15                      # Packing fraction; the particle radius is computed from this value in the unit number density
sigma=0.20                     # Exclusion radius in unit number density; this value must be larger than the particle diameter

# Computational parameters
threads=4                      # Number of threads in openmp
N=200	                         # Number of particles
Nc=5                          # Number of configurations
initial=random                 # Choice of initial condition (random/input)

# Minimization parameters
# - These parameters can be skipped. They are also insensitive in order and case-insensitive.
eps0="tolerance 1e-14"		# energy tolerance for ground states in the unit of v0. Default value is 1e-14.
max_eval="maxsteps 100000" 	# the maximum number of evaluations before quitting calculations. Default value is 10000.
algorithm="algorithm LBFGS"	# Minimization algorithm. Default option is LBFGS. Other available options: LocalGradientDescent, ConjugateGradient, SteepestDescent, MINOP

if [ "$d" == 1 ]; then
v=2.0
elif [ "$d" == 2 ]; then
v=3.1415926535898
elif [ "$d" == 3 ]; then 
v=`echo "4.*${v}/3."| bc -l `
fi 
a=`echo "e(l(${phi_fake}/${v})/${d})" | bc -l`

# Filenames
#exe=${HOME}/code/CCO/CCOptimization/soft_core_stealthy2.out
exe=./MC.out
fname=test

echo "# Run the following command:"
#echo "${exe} ${timelimit} ${eps0} ${max_eval} <<< \"${d} 0 ${Ka} 0.0 1.0 ${sigma} ${phi2} ${threads} ${N} ${Nc} ${Nc_MD} ${initial} ${fname} ground\""

# ground states from Random Initial Conditions
Tfin=1e-20
name=1D_SHU_${chi}_${Tfin}_annealed
# ${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${K1} ${chi} $S0 $vareps0 ${sigma} \
# ${threads} ${N} ${Nc} random $name run " #> log2 #\

${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${K1} ${chi} $S0 $vareps0 ${sigma} \
${threads} ${N} ${Nc} random $name \
stepsize 80.0 adjust 100 equi_steps 10000 \
t_init 0.1 t_fin ${Tfin} coolingrate 0.90 \
acceptance 0.1 0.2 samples_tfin 1 0  run " > ${chi}.log #\
#ground $eps0 $max_eval $algorithm run" > log2

# # ground states from Input Initial Conditions
# loadconfig=${fname}_GS # the name of loaded ConfigPack file.
# K1a=0.				   # change in K1 value.
# ${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${K1a} ${K2a} $S0 $vareps0 ${sigma} ${phi_fake} \
# ${threads} ${N} ${Nc} input ${loadconfig} ${fname}_GS2 \
# ground $eps0 $max_eval $algorithm run" > log4

