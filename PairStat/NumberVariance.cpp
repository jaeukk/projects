/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file NumberVariance.cpp
 *	\brief Function implementations of computing local number variances of point patterns under the periodic boundary conditions. */

#include "NumberVariance.h"
//#include <time.h>


void GetSamplingPoints(DimensionType dim, size_t NumPoints, std::vector<GeometryVector> &centers, RandomGenerator & rng)
{
	centers.clear();
	centers.resize(NumPoints, GeometryVector(static_cast<int> (dim)));

	for (int i = 0; i < NumPoints; i ++){
		for (int j = 0; j < dim; j++){
			centers[i].x[j] = rng.RandomDouble();
		}
	}
}

void GetN(const Configuration &config, std::vector<double> &N, double Rmax, double dR, const GeometryVector &x0_RelativeCoordinates, const std::function<double(const GeometryVector &)> GetRadius)
{
	assert(x0_RelativeCoordinates.Dimension == config.GetDimension()); 
	N.clear();
	N.resize((int)ceil(Rmax / dR), 0.0L);
	//Count the number of particles within a shell of radius (R-dR, R)
	double Rc = sqrt(config.GetDimension())*Rmax;
	config.IterateThroughNeighbors(x0_RelativeCoordinates, Rc, [dR, &N, &GetRadius](const GeometryVector &x0, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void {
		int index = (int)floor(GetRadius(x0) / dR);
		if (index < N.size())
			N[index]++;
	}
	);
	/* cumulative sum */
	for (int i = 1; i < N.size(); i++) {
		N[i] += N[i-1];
	}
}


void MC_NV_single(const Configuration & config, std::vector<GeometryVector> & results, double Rmax, double dR, const std::vector<GeometryVector> & window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct)
{
	results.clear();
	std::vector<double> Ns_0;
	
	config.PrepareIterateThroughNeighbors(Rmax);

	/* First center */
	{
		std::vector<double> N;
		GetN(config, N, Rmax, dR, window_centers[0], WindowSizeFunct);
		results.resize(N.size(), GeometryVector(0,0,0,0));
		Ns_0 = std::vector<double>(N);
		for (int i = 0; i < N.size(); i++){
			results[i].x[0] = dR*(i+1);
		}
	}

	/* second and other centers */
#pragma omp parallel for
	for (int i=1; i< window_centers.size(); i++){
		std::vector<double> N;
		GetN(config, N, Rmax, dR, window_centers[i], WindowSizeFunct);

		for (int j=0; j<N.size(); j++){
			double val = N[j] - Ns_0[j];
#pragma omp atomic
			results[j].x[1] += val;
#pragma omp atomic
			results[j].x[2] += val*val;
		}
	}

	/* summarize data */
	double num_window_centers = (double)window_centers.size();
	for (int i = 0; i < results.size(); i++){
		double mean = results[i].x[1] / num_window_centers;
		double variance = (results[i].x[2]/num_window_centers - mean*mean) ;
		if (num_window_centers > 1.0){
			variance *= num_window_centers/(num_window_centers-1.0) ;
		}

		results[i].x[1] = variance;
		results[i].x[2] = 0.0L;
		results[i].x[3] = 0.0L;
	}
}

void MC_NV_single(const Configuration &config, std::vector<GeometryVector> &results, double Rmax, double dR, size_t num_window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct, int seed)
{
	std::vector<GeometryVector> centers;
	RandomGenerator rng (seed);
	GetSamplingPoints(config.GetDimension(), num_window_centers, centers, rng);
	MC_NV_single(config, results, Rmax, dR, centers, WindowSizeFunct);
}

void MC_NV_Ensemble(const std::function<Configuration(size_t i)> &GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> &Result, double Rmax, double dR, const std::vector<GeometryVector> &WindowCenters, const std::function<double(const GeometryVector &)> WindowSizeFunct)
{
	Result.clear();
	std::cout<< "\nComputing Local Number Variance (given window centers):";
	progress_display pd(NumConfigs);

	/* repeat over configurations */
	for (int i = 0; i < NumConfigs; i++){
		Configuration c = GetConfigsFunction(i);
		std::vector<GeometryVector> Result_;
		MC_NV_single(c, Result_, Rmax, dR, WindowCenters, WindowSizeFunct);
		pd ++ ;
		if (i == 0){
			Result = std::vector<GeometryVector> (Result_);
			for (int j = 0; j < Result.size(); j++){
				Result[j].x[3] = Result[j].x[1];
				Result[j].x[1] = 0.0;
			}
		}
		else{
			for (int j = 0; j < Result.size(); j++){
				double val = Result_[j].x[1] - Result[j].x[3];
				Result[j].x[1] += val;
				Result[j].x[2] += val*val;
			}
		}
	}
	/* Summarize data */
	for(int i = 0; i < Result.size(); i++){
		Result[i].x[1] /= (double)NumConfigs;
		/*SE*/
		if(NumConfigs > 1){
			Result[i].x[2] = std::sqrt(abs(Result[i].x[2]/NumConfigs - Result[i].x[1]*Result[i].x[1])/(NumConfigs - 1));	
		}
		else{
			Result[i].x[2] = 0.0;
		}
		Result[i].x[1] += Result[i].x[3];
		Result[i].x[3] = 0.0;
	}
	std::cout << "Done!\n";
}

void MC_NV_Ensemble(const std::function<Configuration(size_t i)> &GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> &Result, double Rmax, double dR, size_t num_window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct, int seed, bool change_centers)
{
	RandomGenerator rng(seed);
	if(change_centers){
		/* Use different window centers */
		Result.clear();
		std::cout<< "Computing Local Number Variance (varying window centers):";
		progress_display pd(NumConfigs);

		/* repeat over configurations */
		for (int i = 0; i < NumConfigs; i++){
			Configuration c = GetConfigsFunction(i);
			std::vector<GeometryVector> WindowCenters;
			GetSamplingPoints(c.GetDimension(), num_window_centers, WindowCenters, rng);

			std::vector<GeometryVector> Result_;
			MC_NV_single(c, Result_, Rmax, dR, WindowCenters, WindowSizeFunct);
			pd ++ ;
			if (i == 0){
				Result = std::vector<GeometryVector> (Result_);
				for (int j = 0; j < Result.size(); j++){
					Result[j].x[3] = Result[j].x[1];
					Result[j].x[1] = 0.0;
				}
			}
			else{
				for (int j = 0; j < Result.size(); j++){
					double val = Result_[j].x[1] - Result[j].x[3];
					Result[j].x[1] += val;
					Result[j].x[2] += val*val;
				}
			}
		}
		/* Summarize data */
		for(int i = 0; i < Result.size(); i++){
			Result[i].x[1] /= (double)NumConfigs;
			/*SE*/
			if(NumConfigs > 1){
				Result[i].x[2] = std::sqrt(abs(Result[i].x[2]/NumConfigs - Result[i].x[1]*Result[i].x[1])/(NumConfigs-1));	
			}
			else{
				Result[i].x[2] = 0.0;
			}
			Result[i].x[1] += Result[i].x[3];
			Result[i].x[3] = 0.0;
		}
	}
	else{
		/* Use identical window centers */
		std::vector<GeometryVector> window_centers;
		DimensionType d = GetConfigsFunction(0).GetDimension();
		GetSamplingPoints(d, num_window_centers, window_centers,rng);

		MC_NV_Ensemble(GetConfigsFunction, NumConfigs, Result, Rmax, dR, window_centers, WindowSizeFunct);
	}

}

/** \brief A function to generate a random rotation matrix in d dimensions.
*	@param[in] d	Dimension
*	@param[out] R	A random rotation matrix.
*	@param[in] gen	A random number generator. */
// 3D case is incomplete
inline void RandomRotation(DimensionType d, std::vector<GeometryVector> & R, RandomGenerator & gen) {
	R.resize(0);
	double theta = 0;
	double phi = 0;
	double psi = 0;
	switch (d) {
	case 2:
		theta = 2.0*pi*gen.RandomDouble();
		R.push_back(GeometryVector(cos(theta), sin(theta)));
		R.push_back(GeometryVector(-sin(theta), cos(theta)));
		break;
	case 3:
		//***Euler angles***
		phi = 2.0*pi*gen.RandomDouble(); //D
		theta = pi * gen.RandomDouble(); //C
		psi = 2.0*pi*gen.RandomDouble(); //B
										 //R = BCD

		R.push_back(GeometryVector(cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi), cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi), sin(psi)*sin(theta)));
		R.push_back(GeometryVector(-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi), -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi), cos(psi)*sin(theta)));
		R.push_back(GeometryVector(sin(theta)*sin(phi), -sin(theta)*cos(phi), cos(theta)));
		break;
	default:
		std::cout << "NumberVariance::RandomMatrix out of dimension";
		break;
	}
}



void MC_NV_Ensemble_Orientational(const std::function<Configuration(size_t i)> &GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> &Result, double Rmax, double dR, const std::vector<GeometryVector> &WindowCenters, const std::function<double(const GeometryVector &)> WindowSizeFunct)
{
	/* essentially the same as MC_NV_Ensemble(...) function. In each window center, it modifies WindowSizeFunct (...) to account for the random orientation*/
	Result.clear();
	RandomGenerator rng_rot(0);
	DimensionType d = GetConfigsFunction(0).GetDimension();
	std::cout<< "\nComputing Local Number Variance (given window centers):";
	progress_display pd(NumConfigs);

	/* repeat over configurations */
	for (int i = 0; i < NumConfigs; i++){
		Configuration c = GetConfigsFunction(i);
		std::vector<GeometryVector> Result_;
		
		/* Get a rotation matrix */
		std::vector<GeometryVector> R;
		RandomRotation(d, R, rng_rot);
		auto WindowSizeFunct_rotated = [&d, &R, &WindowSizeFunct](const GeometryVector & x)-> double{
			GeometryVector x_ = x;
			for (int i=0; i<d; i++){
				x_.x[i] =0.0;
				for (int j=0; j<d; j++){
					x_.x[i] += R[i].x[j] * x.x[j];
				}
			}
			return WindowSizeFunct(x_);
		};

		MC_NV_single(c, Result_, Rmax, dR, WindowCenters, WindowSizeFunct_rotated);
		pd ++ ;
		if (i == 0){
			Result = std::vector<GeometryVector> (Result_);
			for (int j = 0; j < Result.size(); j++){
				Result[j].x[3] = Result[j].x[1];
				Result[j].x[1] = 0.0;
			}
		}
		else{
			for (int j = 0; j < Result.size(); j++){
				double val = Result_[j].x[1] - Result[j].x[3];
				Result[j].x[1] += val;
				Result[j].x[2] += val*val;
			}
		}
	}
	/* Summarize data */
	for(int i = 0; i < Result.size(); i++){
		Result[i].x[1] /= (double)NumConfigs;
		/*SE*/
		if(NumConfigs > 1){
			Result[i].x[2] = std::sqrt(abs(Result[i].x[2]/NumConfigs - Result[i].x[1]*Result[i].x[1])/(NumConfigs - 1));	
		}
		else{
			Result[i].x[2] = 0.0;
		}
		Result[i].x[1] += Result[i].x[3];
		Result[i].x[3] = 0.0;
	}
	std::cout << "Done!\n";

}
