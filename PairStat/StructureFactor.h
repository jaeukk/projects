/**
 *	Author	: Ge Zhang
 *	Email	: 
 *	Date	:  */



#ifndef STRUCTUREFACTOR_INCLUDED
#define STRUCTUREFACTOR_INCLUDED

#include <functional>
#include <vector>

/* header fileds in $(cores)*/
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <etc.h>
/** \file StructureFactor.h
 *	\brief Header file for computing the structure factor of point configurations,
	spectral density of sphere packings, and spectral density for packings of nonspherical particles (optional). */





//calculate the Structure factor of certain k
double StructureFactor(const Configuration & Config, const GeometryVector & k);
//calculate the Structure factof with weights; here config.getCharacteristics(i) = the weight on the i-th particle
double WeightedStructureFactor(const PeriodicCellList<std::complex<double>> & Config, const GeometryVector & k);
//calculate structure factor of all ks in the range CircularKMax (i.e. for all ks that abs(k)<CircularKMax) and all of their multiples in the range LinearKMax
//Results has the following form:
//( abs(k), S(k), KPrecision, \delta S(k) )
//GetConfigsFunction should return configurations of index i when called
//if SampleProbability is <1, then for any k point satisfying the above condition there is SampleProbability probability that it will be used to calculate S(k) and (1-SampleProbability) probability that it will not be used.
void IsotropicStructureFactor(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision=0.01, double SampleProbability=1.0, size_t option = 0, double CircularKMin = 0.0);

void IsotropicStructureFactor_weighted(std::function<const PeriodicCellList<std::complex<double>>(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision = 0.01, double SampleProbability = 1.0);
//calculate structure factor of all ks in the cubic box [-limit, limit]^d.
//dk : structure factors at ks in [k,k+dk]^d will be averaged: if dk = 1, structure factors are not averaged
//return coordinates of the Grid
//For contourplot in Matlab
std::vector<GeometryVector> DirectionalStructureFactor(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, long limit, std::vector<GeometryVector> & Results, int dk=1, double SampleProbability = 1.0);



#endif