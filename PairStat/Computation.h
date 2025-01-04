/**
 *	Editor	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: March 2022 */

/** \file Computation.h 
 * \brief Header file for a Computation class.
 * It simplified the input and output for some calculations involving an ensemble of point configurations.
 * Computation PairStatisticsCLI code written by Ge.*/


#ifndef COMPUTATION_H__
#define COMPUTATION_H__

#include <set>
#include <cmath>
#include <list>
#include <vector>
#include <omp.h>

/* header files in $(cores) */
#include <etc.h>
#include <PeriodicCellList.h>
#include <GeometryVector.h>

#include "PairCorrelation.h"
#include "StructureFactor.h"
#include "NumberVariance.h"
#include "ContactStatistics.h"

class Computation
{
public:
	int num_threads = 1;

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig) =0;
	virtual void Write(const std::string OutputPrefix) =0;
	virtual void Plot(const std::string OutputPrefix, const std::string & Title) =0;

	//some computations allow additional options
	//any command not understood by PairStatisticsCLI() is passed to here
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output)
	{
		//by default, accept no additional option.
		output << "Unrecognized command!\n";
	}

	virtual ~Computation()
	{
	}

	void SetNumThreads(int n){
		this->num_threads = n;
	}
};

/* class to compute pair correlation functions as a function of radius $r$ */
class g2Computation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double g2rmax;
	double resolution;
	HistogramGenerator * pHGen;
	
	g2Computation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"g2 R_max=";
		ifile>>g2rmax;
		resolution=1.0;
		pHGen = nullptr;
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	virtual ~g2Computation()
	{
		if (pHGen != nullptr)
			delete pHGen;
	}
};

/* class to compute pair correlation functions as a function of a displacement vector */
class Directionalg2Computation : public Computation
{
public:
	double g2rmax;
	size_t NumThetaBins;
	std::vector< std::vector<GeometryVector> > results;

	Directionalg2Computation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"g2 R_max=";
		ifile>>g2rmax;
		ofile<<"Num Theta Bins=";
		ifile>>NumThetaBins;
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title)
	{
		std::cerr<<"Directionalg2Computation::Plot is not implemented!\n";
	}
};

/* a class to compute the structure factor as a function of a wavenumber */
class SkComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double CircularKMax, LinearKMax, KPrecision;
	double SampleProbability;
	size_t average_option = 0;
	
	SkComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"circular K_max=";
		ifile>>CircularKMax;
		ofile<<"linear K_max=";
		ifile>>LinearKMax;
		ofile<<"K precision (binning width)=";
		ifile>>KPrecision;

		SampleProbability = 1.0;
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};



/* a class to compute the structure factor as a function of a wavenumber */
class NearestNeighbors : public Computation
{
public:
	std::vector<GeometryVector> distribution, min_distribution;
	std::vector<double> minPairs;
	size_t SampleDistanceSize;
	double resolution;
	bool logscale_bin = false;
	bool Output_AllMinPairDistances = false;
	
	NearestNeighbors(std::istream & ifile, std::ostream & ofile)
	{
		ofile << "Sample distance sizes (to generate bins) ";
		ifile >> SampleDistanceSize;
		resolution=1.0;	
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


/* class to compute contact numbers */
class CoordinationNumberComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double rmin, rmax;
	int bins;
	
	CoordinationNumberComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"Coordination numbers between rmin and rmax = ";
		ifile>>rmin;
		ifile>>rmax;
		ofile<<"Number of bins = ";
		ifile>>bins;
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	virtual ~CoordinationNumberComputation()
	{}
};

class CumulativeCoordinationNumberComputation : public Computation
{
public:
	std::vector<GeometryVector> result;
	double rmin, rmax, dr;
	
	CumulativeCoordinationNumberComputation(std::istream & ifile, std::ostream & ofile)
	{
		ofile<<"Cumulative Coordination numbers between rmin and rmax = ";
		ifile>>rmin;
		ifile>>rmax;
		ofile<<"dr = ";
		ifile>>dr;
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);

	virtual ~CumulativeCoordinationNumberComputation()
	{}
};


/* a class to compute the local number variance as a function of a wavenumber */
class LocalNumberVariance : public Computation
{
public:
	std::vector<GeometryVector> result;
	double Rmax, dR;
	size_t num_samp_centers, seed;
	RandomGenerator rng;
	bool use_prt_centers = false;
	LocalNumberVariance(std::istream & ifile, std::ostream & ofile)
	{
		ofile << "Rmax = ";
		ifile >> Rmax;
		ofile << "dR = ";
		ifile >> dR;
		ofile << "number of random sampling centers = ";
		ifile >> num_samp_centers;
		ofile << "random seed = ";
		ifile >> seed;
		rng.seed(seed);
	}

	virtual void Compute(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig);
	virtual void Write(const std::string OutputPrefix);
	virtual void Plot(const std::string OutputPrefix, const std::string & Title);
	virtual void ProcessAdditionalOption(const std::string & option, std::istream & input, std::ostream & output);
};


#endif // COMPUTATION_H__