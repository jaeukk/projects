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


void ContactNumbers(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MinDistance, double MaxDistance, std::vector<GeometryVector> & Result);


#endif