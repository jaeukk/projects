/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	August. 2020 */

 /** \file main.cpp
  *	\brief Function implementations of testing the LocalNumber class. */

#include "../GenerateConfigs.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

#include "../etc.h"
#include "../GeometryVector.h"
#include "../PeriodicCellList.h"
#include "../RandomGenerator.h"
#include "NR_distribution.h"


void Example_Bernoulli(std::istream & ifile, std::ostream & ofile) {
	ofile << "=== Test of the LocalNumber class with Bernoilli distribution === \n";

	double p;	// Expectation value
	ofile << "Expectation value (<1) = ";	ifile >> p;
	if (p > 1.0) {
		std::cerr << "p is too large\n";
		return;
	}

	int rng_seed, num_rn;
	ofile << "Seed for a random number generator = "; ifile >> rng_seed;
	ofile << "The number of data that we will collect = "; ifile >> num_rn;

	// A random number generator follows the Bernoulli distribution 
	RandomGenerator rng(rng_seed);
	auto f = [&p, &rng]()->int {
		return (rng.RandomDouble() < p) ? 1 : 0;
	};

	// collect random values.
	ofile << "\nWe collect random values:\n";
	LocalNumber Bernoulli;
	int X;
	for (int i = 0; i < num_rn; i++) {
		X = f();
		Bernoulli.report(X);
	}
	//Compute mean and its standard error
	double mean = 0, SE_mean;	SE_mean = Bernoulli.GetSE_Mean(mean);
	ofile << "After collecting " << num_rn << " samples, the estimated values of expectation value is \n";
	ofile << "p = " << mean << " +- " << SE_mean << std::endl;
	ofile << "\texpected SE = " << sqrt(p*(1 - p) / Bernoulli.GetTotalCount()) << std::endl;

	ofile << "Please the number of data that we will additionally collect = ";
	ifile >> num_rn;
	for (int i = 0; i < num_rn; i++)
		Bernoulli.report(f());

	SE_mean = Bernoulli.GetSE_Mean(mean);
	ofile << "After collecting " << num_rn << " samples, the estimated values of expectation value is \n";
	ofile << "p = " << mean << " +- " << SE_mean << std::endl;
	ofile << "\texpected SE = " << sqrt(p*(1 - p) / Bernoulli.GetTotalCount()) << std::endl;

	/* Variance	*/
	ofile << "From the collected random numbers, the estimated variance and its standard error is\n";
	double variance = 0, SE_var = 0;
	variance = Bernoulli.GetVariance();
	SE_var = Bernoulli.GetSE_Variance();
	ofile << "Variance = " << variance << " +- " << SE_var << std::endl;
	ofile << "\texpected variance = " << p * (1 - p) << std::endl;
}

void Example_DiscreteNormal(std::istream & ifile, std::ostream & ofile) {
	ofile << "=== Test of the LocalNumber class with discrete normal distribution === \n";

	double sigma, mu;	// Expectation value
	ofile << "Please input parameters for normal distribution \n";
	ofile << "Mean = ";		ifile >> mu;
	ofile << "standard deviation (>0) = ";	ifile >> sigma;

	int rng_seed, num_rn;
	ofile << "Seed for a random number generator = "; ifile >> rng_seed;
	ofile << "The number of data that we will collect = "; ifile >> num_rn;

	//A random number generator follows the normal distribution...
	RandomGenerator rng(rng_seed);
	auto g = [&sigma, &mu, &rng]()->int {
		return (int)round(mu + rng.RandomDouble_normal(sigma)); };

	LocalNumber Normal;
	for (int i = 0; i < num_rn; i++)
		Normal.report(g());

	//Compute statistics
	double mean, SE_mean, variance, SE_var;
	SE_mean = Normal.GetSE_Mean(mean);
	variance = Normal.GetVariance();

	ofile << "After collecting " << Normal.GetTotalCount() << " random values,\n";
	ofile << "mean = " << mean << " +- " << SE_mean << std::endl;
	ofile << "\texpected SE = " << sigma / sqrt(Normal.GetTotalCount()) << std::endl;

	SE_var = Normal.GetSE_Variance();
	ofile << "Bootstrap sample size = 50:\n";
	ofile << "Variance = " << variance << " +- " << SE_var << std::endl;
	ofile << "\texpected variance = " << sigma * sigma << std::endl;

	size_t Num_BootstrapSamples;
	ofile << "Please input larger number of Bootstrap samples";
	ifile >> Num_BootstrapSamples;

	SE_var = Normal.GetSE_Variance(Num_BootstrapSamples);
	ofile << "Bootstrap sample size = " << Num_BootstrapSamples << ":\n";
	ofile << "Variance = " << variance << " +- " << SE_var << std::endl;


	ofile << "\n\n"
		<< "==== File IO ====\n";

	std::string name;
	ofile << "Please input the prefix to store the histogram.\n";
	ifile >> name;
	std::fstream save((name + ".bin").c_str(), std::fstream::out | std::fstream::binary);
	Normal.WriteBinary(save);
	save.close();
	ofile << "Data are saved in " << (name + ".bin").c_str()
		<< "\n\n"
		<< "Now, we open the file to load data in it. \n";
	std::fstream load((name + ".bin").c_str(), std::fstream::in | std::fstream::binary);
	LocalNumber Normal2(load);
	load.close();
}

void Example_NV_Binomial(std::istream & ifile, std::ostream & ofile) {
	ofile << "=== N(R) for binomial point patterns ===\n";

	DimensionType dim = 3;
	size_t NumSampWindow = 1000, NumConfigs = 200;
	double Rmax = 3, dR = 0.01, rho = 27;
	int seed_window, seed_config;

	ofile << "Use MC_nv () and MC_nv_distribution() functions\n"
		<< "Throw " << NumSampWindow << " randomly placed spherical windows of radius " << Rmax << "\n";

	ofile << "Random seed for sampling windows = ";
	ifile >> seed_window;
	ofile << "Random seed for configurations = ";
	ifile >> seed_config;
	RandomGenerator rng(seed_window);
	PoissonPtGen C(dim, 0, seed_config);
	C.SetRho(rho);

	auto GetConfig = [&C](size_t i)->Configuration {
		Configuration c = GetUnitCubicBox(3, 0.1);
		c.Resize(1000);
		C.GenerateC(c);
		return c;
	};

	std::vector<GeometryVector> variances, samp_center;
	samp_center = GetSamplingPoints(3, NumSampWindow, seed_window);

	//Ordinary function to compute the local number variance.
	ofile << "Compute local number variance \n";
	std::string name = "3D_BP/3D_BP_27000";
	ofile << "Data will be save in files with a prefix " << name << std::endl;

	MC_nv(GetConfig, NumConfigs, variances, Rmax, dR, samp_center);
	WriteFunction(variances, (name + "_LNV").c_str());

	//LocalNumber classes
	BootstrapSetup setup(name);
	setup.ComputeVariance = true;
	setup.ComputeSkewness = true;
	setup.ComputeKurtosis = true;
	ofile << "We use bootstrap method to estimate standard error\n";
	MC_nv_distribution(GetConfig, NumConfigs, variances, Rmax, dR, setup, samp_center);
	WriteFunction(variances, (name + "_LNV2").c_str());
	ofile << "end";
}

void Example_NV_Poisson(std::istream & ifile, std::ostream & ofile) {
	ofile << "=== N(R) for Poisson point patterns ===\n";

	DimensionType dim = 3;
	size_t NumSampWindow = 1000, NumConfigs = 200;
	double Rmax = 3, dR = 0.01, rho = 27;
	int seed_window, seed_config;

	ofile << "Use MC_nv () and MC_nv_distribution() functions\n"
		<< "Throw " << NumSampWindow << " randomly placed spherical windows of radius " << Rmax << "\n";

	ofile << "Random seed for sampling windows = ";
	ifile >> seed_window;
	ofile << "Random seed for configurations = ";
	ifile >> seed_config;
	RandomGenerator rng(seed_window);
	PoissonPtGen C(dim, 1, seed_config);
	C.SetRho(rho);

	auto GetConfig = [&C](size_t i)->Configuration {
		Configuration c = GetUnitCubicBox(3, 0.1);
		c.Resize(1000);
		C.GenerateC(c);
		return c;
	};

	std::vector<GeometryVector> variances, samp_center;
	samp_center = GetSamplingPoints(3, NumSampWindow, seed_window);

	//Ordinary function to compute the local number variance.
	ofile << "Compute local number variance \n";
	std::string name = "3D_PPP/3D_PPP_27000";
	ofile << "Data will be save in files with a prefix " << name << std::endl;

	MC_nv(GetConfig, NumConfigs, variances, Rmax, dR, samp_center);
	WriteFunction(variances, (name + "_LNV").c_str());

	//LocalNumber classes
	BootstrapSetup setup(name);
	setup.ComputeVariance = true;
	setup.ComputeSkewness = true;
	setup.ComputeKurtosis = true;
	ofile << "We use bootstrap method to estimate standard error\n";
	MC_nv_distribution(GetConfig, NumConfigs, variances, Rmax, dR, setup, samp_center);
	WriteFunction(variances, (name + "_LNV2").c_str());
	ofile << "end";
}

void Example_NV_CubicLattice(std::istream & ifile, std::ostream & ofile) {
	ofile << "=== N(R) for Cubic lattice ===\n";
	double Rmax = 100, dR = 0.01;
	size_t NumSampWindow_lattice = 1000000;
	int seed_window;
	DimensionType dim = 2;
	ofile << "Use MC_nv () and MC_nv_distribution() functions"
		<< "to investigate distribution of N(R) for cubic lattice\n"
		<< "Throw " << NumSampWindow_lattice << " randomly placed spherical windows of radius " << Rmax << "\n";

	ofile << "Random seed for sampling windows = ";
	ifile >> seed_window;
	//Define a configuration for Cubic lattice.
	LatticeGen Latt("Square");
	Latt.SetNumCells(1);
	auto GetConfig = [&Latt](size_t i)->Configuration {
		Configuration c;
		Latt.GenerateC(c);
		return c;
	};
	std::vector<GeometryVector> variances;
	std::vector<GeometryVector> samp_center = GetSamplingPoints(dim, NumSampWindow_lattice, seed_window);

	ofile << "Compute local number variance \n";
	//std::string name = "3D_Cubic/3D_Cubic";
	std::string name = "2D_Square/2D_Square";
	ofile << "Data will be save in files with a prefix " << name << std::endl;

	MC_nv(GetConfig, 1, variances, Rmax, dR, samp_center);
	WriteFunction(variances, (name + "_LNV").c_str());

	//LocalNumber classes
	BootstrapSetup setup(name);
	setup.ComputeVariance = true;
	setup.ComputeSkewness = true;
	setup.ComputeKurtosis = true;
	ofile << "We use bootstrap method to estimate standard error\n";
	MC_nv_distribution(GetConfig, 1, variances, Rmax, dR, setup, samp_center);
	WriteFunction(variances, (name + "_LNV2").c_str());
	ofile << "end";
}

void Example_NV_RSA(std::istream & ifile, std::ostream & ofile) {
	ofile << "=== N(R) for saturated RSA packings ===\n";

	DimensionType dim = 3;
	size_t NumSampWindow = 1000, NumConfigs = 200;
	double Rmax = 10, dR = 0.01;
	int seed_window, seed_config;

	ofile << "Use MC_nv () and MC_nv_distribution() functions\n"
		<< "Throw " << NumSampWindow << " randomly placed spherical windows of radius " << Rmax << "\n";

	ofile << "Random seed for sampling windows = ";
	ifile >> seed_window;
	ofile << "Random seed for configurations = ";
	ifile >> seed_config;

	RSA_Gen C(dim, 4 * 4 * 4 * 1000, 1.0);
	C.SetNumThreads(omp_get_num_threads());
	auto GetConfig = [&C](size_t i)->Configuration {
		Configuration c;
		c.Resize(1000);
		C.GenerateC(c);
		return c;
	};

	std::vector<GeometryVector> variances, samp_center;
	samp_center = GetSamplingPoints(3, NumSampWindow, seed_window);

	//Ordinary function to compute the local number variance.
	ofile << "Compute local number variance \n";
	std::string name = "3D_RSA/3D_RSA_sat_64000";
	ofile << "Data will be save in files with a prefix " << name << std::endl;

	MC_nv(GetConfig, NumConfigs, variances, Rmax, dR, samp_center);
	WriteFunction(variances, (name + "_LNV").c_str());

	//LocalNumber classes
	BootstrapSetup setup(name);
	setup.ComputeVariance = true;
	setup.ComputeSkewness = true;
	setup.ComputeKurtosis = true;
	ofile << "We use bootstrap method to estimate standard error\n";
	MC_nv_distribution(GetConfig, NumConfigs, variances, Rmax, dR, setup, samp_center);
	WriteFunction(variances, (name + "_LNV2").c_str());
	ofile << "end";
}
