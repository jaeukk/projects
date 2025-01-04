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

void CumulativeCoordinationNumber(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double r_min, double r_max, double dr, std::vector<GeometryVector> &Result)
{
//	double log_ratio = log(MaxDistance/MinDistance)/(double)num_bins;

	/* Determine the min distances*/
	double rmin = 100000000.0;
	#pragma omp parallel for 
	for (size_t i = 0; i < NumConfigs; i++){
		Configuration c = GetConfigsFunction(i);
		double temp = 100000000.0;
		for(size_t j = 0; j < c.NumParticle(); j++){
			temp = std::min(temp, c.NearestParticleDistance(j));
		}
		#pragma omp critical
		{
			rmin = std::min(rmin, temp);
		}
	}
	if (rmin > r_min){
		std::cout << "r_min parameter is smaller than the minimum distance, so we take r_min to be the shortest distance\n";
	}
	else{
		rmin = r_min;
	}

	/* Initialize the Result vector. */
	Result.clear();
	size_t num = (size_t) ceil((r_max-rmin)/dr);
	{
		double r = rmin;
		for (size_t i=0; i<num; i++){
			Result.emplace_back(r,0,0);
			r+= dr;
		}
	}
	

	if(Verbosity>2)
		std::cout<<"computing Cumulative Coordination Number Z(r)";
	progress_display pd(NumConfigs);

	for(size_t i=0; i<NumConfigs; i++){
		Configuration c = GetConfigsFunction(i);
		c.PrepareIterateThroughNeighbors(r_max);
		std::vector<double> num_list_config(Result.size(),0.0);

		#pragma omp parallel for 
		for (size_t j=0; j<c.NumParticle(); j++)
		{
			std::vector<double> num_list(Result.size(),0.0);
			c.IterateThroughNeighbors(
				c.GetRelativeCoordinates(j), r_max, 
				[&num_list, j, r_max, num, dr, rmin](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void{
					double distance = sqrt(shift.Modulus2());
					size_t idx = 0;
					if (SourceAtom != j){
						if (distance > rmin){
							idx = (size_t) ceil((distance-rmin)/dr);
						}
						if (idx < num_list.size()){
							num_list.at(idx) ++;
						}
					}
				}
			);
			for (size_t k=1; k<num_list.size(); k++){
				num_list[k] += num_list[k-1];
			}
			
			for (size_t k=0; k < num_list.size(); k++){
				#pragma atomic
				num_list_config[k] += num_list[k];
			}
		}

		for (size_t k=0; k < num_list_config.size(); k++){
			double val = num_list_config[k]/c.NumParticle();
			Result[k].x[1] += val;
			Result[k].x[2] += val *val ;
		}

		pd ++ ;
	}

	/* Summary */
	for (size_t k=0; k < Result.size(); k++){
		Result[k].x[1] /= NumConfigs;
		double var = (Result[k].x[2]/NumConfigs - Result[k].x[1]*Result[k].x[1])/(NumConfigs-1.0);
		Result[k].x[2] = (var > 0.0)? sqrt(var) : 0.0;
	}
	
	if(Verbosity>2)
		std::cout<<"Done\n";
}
