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

double CumulativeCoordinationNumber(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double r_min, double r_max, double dr, std::vector<GeometryVector> &Result, bool logscale)
{
//	double log_ratio = log(MaxDistance/MinDistance)/(double)num_bins;

	/* Determine the min distance*/
	double rmin = 100000000.0;
	double ratio = 1.01;
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
		std::cout << "Parameter r_min is smaller than the actual minimum distance. We now take r_min to be the shortest distance (" << rmin << ")\n";
	}
	else{
		rmin = r_min;
	}

	/* Initialize the Result vector. */
	Result.clear();
	std::function<size_t(double)> get_bin = nullptr;

	if (logscale){
		get_bin = [ratio,rmin,dr](double r) -> size_t {
			if (r <= rmin){
				return 0;
			}
			else{
				return (size_t)std::ceil(std::log(1 + (ratio - 1.)*(r-rmin)/dr) / std::log(ratio)) - 1;
			}
		};
		std::cout << "Generate log-scaled bins.\n";
		std::cout << "bin sizes follow a geometric sequence: delta r_i = " << dr << " * " << ratio << " ^ i \n";
		double bin = dr;
		double r = rmin + dr; 
		size_t num = get_bin (r_max); // (size_t) ceil( std::log( (ratio-1.0)*(r_max-rmin)/dr + 1.) / std::log(ratio));
		Result.reserve(num);
		for (; r < r_max; r += bin){
			Result.emplace_back(r, 0, 0);
			bin *= ratio;
		}
		std::cout << "There are " << num << " bins \n";
	}
	else{
		get_bin = [rmin,dr](double r) -> size_t {
			if (r <= rmin){
				return 0;
			}
			else{
				return (size_t)std::floor((r-rmin)/dr);
			}
		};
		size_t num = get_bin(r_max);
		{
			double r = rmin;
			for (size_t i=0; i<num; i++){
				r+= dr;
				Result.emplace_back(r,0,0);
			}
		}
	}
	

	if(Verbosity>2)
		std::cout<<"computing Cumulative Coordination Number Z(r)";
	progress_display pd(NumConfigs);

	for(size_t i = 0; i < NumConfigs; i++){
		Configuration c = GetConfigsFunction(i);
		c.PrepareIterateThroughNeighbors(r_max);
		std::vector<double> num_list_config(Result.size(),0.0);  // for configuration i

		#pragma omp parallel for shared(num_list_config)
		for (size_t j = 0; j < c.NumParticle(); j++)
		{
			std::vector<double> num_pairs(Result.size(),0.0); // from particle j
			c.IterateThroughNeighbors(
				c.GetRelativeCoordinates(j), r_max, 
				[&num_pairs, rmin, j, &get_bin](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void{
					if (SourceAtom != j){
						double distance = sqrt(shift.Modulus2());
						size_t idx = get_bin(distance);
						if (idx < num_pairs.size()){
							num_pairs.at(idx) ++;
						}
					}
				}
			);
			/* cumulative sum */
			// for (size_t k = 1; k < num_list.size(); k++){
			// 	#pragma atomic
			// 	num_list.at(k) += num_list[k-1];
			// }
			
			for (size_t k = 0; k < num_pairs.size(); k++){
				#pragma omp atomic
				num_list_config.at(k) += num_pairs[k];	// now, num_list_config =  # of all pairs of a certain distance in a configuration.
			}
		}

		double CumSum = 0.0L, val = 0.0L;
		for (size_t k = 0; k < num_list_config.size(); k++){
			
			// if ( k == 0 ){
				// 	val = num_list_config[k] / c.NumParticle();
				// }
				// else{
				// 	/* cumulative sum */
				// 	val += num_list_config[k] / c.NumParticle();
				// }
					
			CumSum += num_list_config[k];
			val = CumSum / c.NumParticle();
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
	return rmin;
}
