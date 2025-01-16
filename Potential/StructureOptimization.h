/** author: Ge Zhnag
	date:
	email:		
	comments: Jaeuk Kim (June 2019).	*/

/**	\file StructureOptimization.h
	\brief Header file to perform energy minimization... */


#ifndef STRUCTUREOPTIMIZATION_INCLUDED
#define STRUCTUREOPTIMIZATION_INCLUDED

#include <cassert>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <nlopt.h>
/* in $(cores) */
#include <PeriodicCellList.h>

#include "Potential.h"

//default: 0
//if not zero, then partially optimized configurations are saved every these steps
//will save to two files: "StructureOptimization_Progress_a.configuration" 
//and "StructureOptimization_Progress_b.configuration" alternate
extern size_t StructureOptimization_SaveConfigurationInterval;
extern bool StructureOptimization_SaveNumberedConfiguration;

class NotANumberFound : public std::exception
{
public:
	NotANumberFound()
	{
	}
	virtual const char * what() const throw()
	{
		return "Caught a NaN!\n";
	}
};

/** \brief Minimize the potential energy of a given initial configuration via MINOP algorithm. 
 *	@param[in,out]	List	An initial configuration, and it becomes optimized.
 *	@param[in]		pot		A refernece of the potential energy.
 *	@param[in]		Pressure	A constant pressure in the isobaric ensembles.
 *	@param[in]		Switch	Determine which one is being optimized;
					0:	particle positions.
					1:	fundamental cell
					2:	both particle positions and fundamental cell
 *	@param[in]		MinDistance	Particle displacement is smaller than this value, the movement is rejected.
 *	@param[in]		MaxStep	The maximal number of optimization steps.
 *	@param[in]		Emin	*/
void RelaxStructure_MINOP_withoutPertubation(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000, double Emin=-1e10);

//do MINOP relax without pertubation, then with pertubation, return the best one
void RelaxStructure_MINOP(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000);

//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_NLOPT(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_NLOPT_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double Emin);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_SteepestDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
void RelaxStructure_ConjugateGradient(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000, double minG=1e-13);
void RelaxStructure_LocalGradientDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);

//slower but does not need derivative information
void RelaxStructure_NLOPT_NoDerivative(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, double rescale=1.0, size_t MaxStep=100000, double fmin=(-1.0)*MaxEnergy);

inline void RelaxStructure(Configuration & List, Potential & pot, double Pressure, double MinDistance)
{
	try
	{
		RelaxStructure_NLOPT(List, pot, Pressure, 2, MinDistance);
	}
	catch(EnergyDerivativeToBasisVectors_NotImplemented a)
	{
		for(int i=0; i<5; i++)
		{
			RelaxStructure_NLOPT(List, pot, Pressure, 0, MinDistance);
			RelaxStructure_NLOPT_NoDerivative(List, pot, Pressure, 1, MinDistance);
		}
		RelaxStructure_NLOPT_NoDerivative(List, pot, Pressure, 2, MinDistance);
	}
}

//calculate Elastic Constants
//relax the structure before calling this function
//\epsilon_{ij}=C_{ijkl}*e_{kl}
double ElasticConstant(Potential & pot, Configuration structure, DimensionType i, DimensionType j, DimensionType k, DimensionType l, double Pressure, bool InfiniteTime=false);

//if presult is not nullptr, then presult[0..dimension^4] will be filled with C_{ijkl} 
void PrintAllElasticConstants(std::ostream & out, Configuration stru, Potential & pot, double Pressure, bool InfiniteTime = false, double * presult=nullptr);


//optimize the parameters for a elastic model and rotations so that it fits the target
void ElasticOptimize(std::function<BigVector(const std::vector<double> & param)> getElasticFunc, int nparam, std::vector<double> & param, std::vector<double> ubound, std::vector<double> lbound, const Configuration & c, Potential & pot, double pressure);
#endif
