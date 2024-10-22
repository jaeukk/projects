/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: June 2024 */


/** \file ContactStatistics.h 
 * \brief Header file for a 
 * */

#include "ContactStatistics.h"


void ContactNumbers(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, 
	double MinDistance, double MaxDistance, std::vector<GeometryVector> & Result){

	if(!(MaxDistance>0.0))
		return;
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"computing Local Coordination Numbers";
	progress_display pd(NumConfigs);

	Result.clear();
	{
		Configuration c = GetConfigsFunction(0);
		DimensionType d = c.GetDimension();
		for (int i=0; i<4*d+2; i++){
			Result.emplace_back(i, 0.,0.,0.);
		}
	}

	double mean = 0.0;
	double Min2 = MinDistance*MinDistance;
	#pragma omp parallel for reduction(+:mean)
	for (size_t i=0; i<NumConfigs; i++){
		Configuration c = GetConfigsFunction(i);
		c.PrepareIterateThroughNeighbors(MaxDistance);
		double increment = 1.0/ c.NumParticle();
		std::vector<double> Neighbors (Result.size(), 0);
		double mean_config = 0.0;
		for (size_t j=0; j < c.NumParticle(); j++){
			size_t count = 0;
			c.IterateThroughNeighbors(
				c.GetRelativeCoordinates(j), MaxDistance, 
					[&count, &j, Min2](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void{
						if ( (SourceAtom != j) && (shift.Modulus2() > Min2)){
							count ++ ;
						}
					}
			);
			if (count < Result.size()){
				Neighbors[count] += increment;
			}
			else{
				Neighbors[Result.size()-1] += increment;
			}
			mean_config += (double)count;			
		}
		mean_config /= c.NumParticle();

		mean += mean_config;

		for (size_t j = 0; j<Result.size(); j++){
		#pragma omp atomic 
			Result[j].x[1] += Neighbors[j];
		#pragma omp atomic
			Result[j].x[2] += Neighbors[j]*Neighbors[j];
		}

		#pragma omp critical
			pd ++ ;
		
	}

	mean /= (double)NumConfigs;
	for (size_t j = 0; j<Result.size(); j++){
		double val = Result[j].x[1] / NumConfigs;
		Result[j].x[1] = val;
		Result[j].x[2] = std::sqrt( std::abs(Result[j].x[2] / NumConfigs - val*val) / (NumConfigs-1));
	}

	if(Verbosity>2)
		std::cout<<"Mean Coordination Number = " << mean << std::endl;
}
