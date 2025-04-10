/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: June 2024 */


/** \file ContactStatistics.h 
 * \brief Header file for a 
 * */

#ifndef Z_INCLUDED
#define Z_INCLUDED

#include "PeriodicCellList.h"
#include <vector>
#include <omp.h>
#include <math.h>


/// @brief Measure a distribution of coordination numbers in a range between MinDistance and MaxDistance.
/// @param GetConfigsFunction 
/// @param NumConfigs 
/// @param MinDistance 
/// @param MaxDistance 
/// @param Result 
void ContactNumbers(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MinDistance, double MaxDistance, std::vector<GeometryVector> & Result);

/// @brief Measure cumulative coordination number as a function of distance.
/// @param GetConfigsFunction 
/// @param NumConfigs 
/// @param r_min 
/// @param r_max 
/// @param dr 
/// @param Result 
double CumulativeCoordinationNumber(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double r_min, double r_max, double dr, std::vector<GeometryVector> &Result, bool logscale=false);

#endif