/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	August. 2020 */

 /** \file main.cpp
  *	\brief Function implementations of computing the distribution of N(R)
  * from a ConfigPack files. */



#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

// In $core
#include <etc.h>
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>
size_t Verbosity = 2;


#include "NR_distribution.h"

void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields);
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields);
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns);
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns);

/* Please compile the "examples.cpp" file to use the following functions. */
//Analyze the empirical distribution from random numbers follow the Bernoulli distribution
//of expectation value p.
void Example_Bernoulli(std::istream & ifile, std::ostream & ofile);

//Analyze the empirical distribution from random numbers follow the discrete normal distribution
//of expectation value p.
void Example_DiscreteNormal(std::istream & ifile, std::ostream & ofile);

void Example_NV_Binomial(std::istream & ifile, std::ostream & ofile);
void Example_NV_Poisson(std::istream & ifile, std::ostream & ofile);
void Example_NV_CubicLattice(std::istream & ifile, std::ostream & ofile);
void Example_NV_RSA(std::istream & ifile, std::ostream & ofile);

//#include "../GenerateConfigs.h"

/*	Analyze an ConfigurationPack with the LocalNumber class,
 * or load the data. */
int main(int argc, char ** argv) {
	std::string tempStr;
	char tempstring[1000];
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	RandomGenerator rng(0);
	omp_set_num_threads(4);

	// {
	// 	ofile << "----- 1D saturated RSA ----\n";
	// 	int num = 100, numc=100, numthreads=4;
	// 	ofile << "number of particles = ";
	// 	ifile >> num;
	// 	ofile << "\nnumber of configurations = ";
	// 	ifile >> numc;

	// 	RSA_Gen rsa(1, num, 0.75, 1);
	// 	rsa.SetNumThreads(numthreads);

	// 	for (int i=0; i<numc; i++){
	// 		SpherePacking c;
	// 		rsa.GenerateP(c);
	// 		sprintf(tempstring, "./1D-rsa-0.75_N-%d_%d", num, i);
	// 		WriteConfiguration(c, tempstring);

	// 	}

	// }
	// return 0;


	{
		double dx = -0.1;
		bool printPDF = false, printCDF = false, periodic_boundary = true;
		if (argc >= 3) {
			printPDF = (std::atoi(argv[1]) == 1);
			printCDF = (std::atoi(argv[2]) == 1);
		}
		if (argc >= 4) {
			periodic_boundary = (std::atoi(argv[3]) == 1);
		}
		if (argc == 5){
			dx = std::atof(argv[4]);
		}
		//printPDF = true;
		DimensionType dim = 3;
		size_t NumSampWindow = 1000, startIndex = 0, NumConfigs = 200;
		int seed_window, num_threads;

		double Rmax = 3, dR = 0.01, rho = 27;

		/* set configurations */
		ofile << "Configurations to analyze = ";
		ifile >> name;
		ofile << "\nstart index = ";
		ifile >> startIndex;
		ofile << "\nthe number of configurations = ";
		ifile >> NumConfigs;

		ConfigurationPack configs(name);
		if (configs.NumConfig() < NumConfigs + startIndex) {
			NumConfigs = configs.NumConfig() - startIndex;
			ofile << "The number of available configurations is " << configs.NumConfig() << " in total.\n"
				<< "NumConfigs is set to be " << NumConfigs << "\n";
		}

		auto GetConfig = [&configs, &startIndex](size_t i)-> Configuration {
			return configs.GetConfig(startIndex + i);
		};

		dim = GetConfig(0).GetDimension();


		/* settings for sampling windows */
		ofile << "\n\n"
			<< "-------------------------------------\n"
			<< "\t Settings for sampling windows\n"
			<< "-------------------------------------\n";
		ofile << "Random seed for sampling windows = ";
		ifile >> seed_window;

		ofile << "\nThe number of sampling windwos = ";
		ifile >> NumSampWindow;
		ofile << "\nMaximum radius = ";
		ifile >> Rmax;
		ofile << "\nchange in radius = ";
		ifile >> dR;

		/* Initialize vectors. */
		std::vector<GeometryVector> variances, /* local number variance as a function of R. */
			samp_center;	/* centers on sampling windows in relative coordinates. */

		if (periodic_boundary) {
			ofile << "** Periodic boundary condition is on\n";
			if (dx < 0.0){
				samp_center = GetSamplingPoints(dim, NumSampWindow, seed_window);
			}
			else{
				/* renew the number of sampling windows*/
				int Nx = (int)std::round(std::pow(NumSampWindow, 1.0/dim));
				NumSampWindow = (int) std::pow(Nx,dim);
				samp_center.resize(NumSampWindow, GeometryVector(static_cast<int> (dim)));

				std::vector<double> L;
				{
					Configuration c_temp = GetConfig(0);
					for (int i=0; i<dim; i++)
						L.emplace_back(std::sqrt(c_temp.GetBasisVector(i).Modulus2()));
				}

#pragma omp parallel for 
				for (int i=0; i<samp_center.size(); i++){
					int val = i, index=0;
					for (int j = 0; j<dim; j++){
						index = val % Nx;
						val /= Nx;
						samp_center[i].x[j] = dx*(double)index / L[j];
					}
				}
			}
		}
		else {
			ofile << "** Periodic boundary condition is off (square or cubic boxs only)\n";
			ofile << "** Side length of a box should be larger than 2*R_max\n";
			/* Window centers (samp_center) must lie in a closed cube [R_max, L-R_max)^d. 
			 * To do so, we carry out the following transform:
			 *		x_i -> (1.0-2.0*Rmax/L)* (x_i -0.5),
			 * where x is a vector in samp_center, x_i is the i-th component of x,
			 * and L is the box size.	 */
			Configuration c_temp = GetConfig(0);
			if (dx < 0.0){
				samp_center = GetSamplingPoints(dim, NumSampWindow, seed_window);

				GeometryVector x_center(static_cast<int>(c_temp.GetDimension()));
				for (int i = 0; i < x_center.Dimension; i++) {
					x_center.x[i] = 0.5;
				}
				double factor = 1.0 - 2.0*Rmax / std::sqrt(c_temp.GetBasisVector(0).Modulus2());
				if (factor < 0)
					std::cerr << "2*Rmax must be smaller than simulation box\n";

	#pragma omp parallel for 
				for (int i = 0; i < samp_center.size(); i++) {
					samp_center[i] = samp_center[i] - x_center;	
					samp_center[i] = factor * samp_center[i];
					samp_center[i] = samp_center[i] + x_center;
				}
			}
			else{
				/* renew the number of sampling windows*/
				int Nx = (int)std::round(std::pow(NumSampWindow, 1.0/dim));
				NumSampWindow = (int) std::pow(Nx,dim);
				samp_center.resize(NumSampWindow, GeometryVector(static_cast<int> (dim)));

				std::vector<double> L;
				{
					Configuration c_temp = GetConfig(0);
					for (int i=0; i<dim; i++)
						L.emplace_back(std::sqrt(c_temp.GetBasisVector(i).Modulus2()));
				}

#pragma omp parallel for 
				for (int i=0; i<samp_center.size(); i++){
					int val = i, index=0;
					for (int j = 0; j<dim; j++){
						index = val % Nx;
						val /= Nx;
						samp_center[i].x[j] = (Rmax+dx*(double)index) / L[j];
					}
				}
			}
		}
		ofile << "\nthread number = ";
		ifile >> num_threads;
		omp_set_num_threads(num_threads);

		ofile << "\n\n"
			<< "name of data storage = ";
		ifile >> name;
		//LocalNumber classes
		BootstrapSetup setup(name);
		setup.BootstrapSize = 5;
		setup.ComputeVariance = true;
		setup.ComputeSkewness = true;
		setup.ComputeKurtosis = true;
		setup.ComputeKL = true;//true;
		setup.ComputeL2norm = true;//true;

		setup.SaveData = true;
		setup.PrintPDF = printPDF;
		setup.PrintCDF = printCDF;
		//ofile << "We use bootstrap method to estimate standard error\n";

		std::string save_prefix;
		ofile << "\nprefix for the output files = ";
		ifile >> save_prefix;
		ofile << "\n\n";
		setup.SetOutputPrefix(save_prefix);

		ofile << "Print PDF of distribution = " << setup.PrintPDF << "\n";
		ofile << "Print CDF of distribution = " << setup.PrintCDF << "\n";

		//MC_nv_distribution_PostSampling(GetConfig, NumConfigs, variances, Rmax, dR, setup, samp_center);
		MC_nv_distribution_PostSampling(GetConfig, NumConfigs, variances, Rmax, dR, setup, samp_center);
		ofile << "end";


	}

	return 0;

}


