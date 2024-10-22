/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	August. 2020 */

 /** \file NR_distribution.h
  *	\brief Function implementations of computing the distribution of N(R). */


#include "NR_distribution.h"

  //--------------------------------------------
  /** Member functions of LocalNumber class. */
  //--------------------------------------------


  /** \brief Give an index of random number according to the empirical distribution.*/
inline size_t LocalNumber::GetRandomIndex(RandomGenerator & rng) const {
	double probability = rng.RandomDouble();
	// std::lower_bound ( ..., probability ) gives an iterator of the smallest element that is not smaller than probability.
	return std::lower_bound(this->CDF.begin(), this->CDF.end(), probability) - this->CDF.begin();
}

/** A member function to compare all counts.*/
bool LocalNumber::CheckSum() {
	size_t sum = 0;
	for (auto c = counts.begin(); c != counts.end(); c++)
		sum += *c;

	return sum == totalCount;
}


/** Get an estimated mean of this distribution. */
double LocalNumber::GetMean() const {
	double mean = 0;// guessed_av = 0.5*(this->min + this->max);
#pragma omp parallel for reduction(+:mean)
	for (int i = 0; i < this->counts.size(); i++)
		mean += ((double)i) * (double)this->counts[i] / this->totalCount;

	return this->min + mean;
}
/** Get a standard error in the estimated mean. */
double LocalNumber::GetSE_Mean(double & mean) const {
	mean = this->GetMean();
	double SE = 0.0, mean_i = mean - this->min;
#pragma omp parallel for reduction(+:SE)
	for (int i = 0; i < this->counts.size(); i++)
		SE += i * i*(double)this->counts[i] / this->totalCount;

	SE = (SE - mean_i * mean_i) / (this->totalCount - 1);
	SE = (SE > 0) ? sqrt(SE) : 0.0;
	return SE;
}
/** Get an estimated variance. */
double LocalNumber::GetVariance() const {
	double mean = this->GetMean() - this->min;	//Mean of the indices.
	double variance = 0.0;
#pragma omp parallel for reduction(+:variance)
	for (int i = 0; i < this->counts.size(); i++) {
		double dx = (double)i - mean;
		variance += dx * dx*(double)this->counts[i] / this->totalCount;
	}

	variance *= (1.0 + 1.0 / (this->totalCount - 1));
//	variance = (variance - mean * mean) * (1.0 + 1.0 / (this->totalCount - 1));
	return std::max(0.0, variance);
}
/** Get a standard error in the estimated variance.
 * Bootstrap method is used.
 * @param Bootstrap_Size	The number of bootstrap replica.
 * @return Standard error of variance via the Bootstrap algorithm. */
double LocalNumber::GetSE_Variance(int Bootstrap_Size, int seed) const {
	double variance = 0.0, SE = 0.0;
	UpdateCDF();

#pragma omp parallel for schedule (guided)
	for (int i = 0; i < Bootstrap_Size; i++) {
		std::unique_ptr<LocalNumber> sample = this->GetBootstrapSample(i);
		double temp = sample->GetVariance();
#pragma omp atomic
		variance += temp;
#pragma omp atomic
		SE += temp * temp;
	}

	variance /= Bootstrap_Size;
	SE = (SE / Bootstrap_Size - variance * variance)*(1.0 + 1.0 / (Bootstrap_Size - 1));	//This qunatity sometimes can be negative due to numerical precision...
	SE = (SE > 0) ? sqrt(SE) : 0.0;
	return SE;
}

/** Compute CDF again if CDF_uptodate = false. */
void LocalNumber::UpdateCDF() {
	// CDF is not updated yet. Recompute CDF from this->counts.
	if (!CDF_uptodate) {
		// this->counts is not empty
		if (counts.size() > 0) {
			// The number of bins changed (it always increases) => increases the size of CDF.
			if (max - min + 1 > CDF.size()) {
				this->CDF.resize(max - min + 1);
			}
			//Recompute this->CDF
			size_t i = 0;
			CDF[0] = (double)counts[0] / totalCount;
			for (i = 1; i < counts.size() - 1; i++) {
				CDF[i] = CDF[i - 1] + (double)counts[i] / totalCount;
			}
			//The final value of CDF is always 1.
			CDF[i] = 1.0;
			CDF_uptodate = true;
		}
		else {
			std::cerr << "This object is empty.\n";
		}
	}
}

void LocalNumber::UpdateCDF() const {
	if (!CDF_uptodate) {
		if (counts.size() > 0) {
			if (max - min + 1 > CDF.size()) {
				this->CDF.resize(max - min + 1);
			}
			size_t i = 0;
			CDF[0] = (double)counts[0] / totalCount;
			for (i = 1; i < counts.size() - 1; i++) {
				CDF[i] = CDF[i - 1] + (double)counts[i] / totalCount;
			}
			CDF[i] = 1.0;
			CDF_uptodate = true;
		}
		else {
			std::cerr << "This object is empty.\n";
		}
	}
}



/** Give a histogram.
 * @param result	A vector of GeometryVectors that describe the table of histogram.
		result[i].x[0] = x value
		result[i].x[1] = probability density
		result[i].x[2] = error in x
		result[i].x[3] = error in y (via multinomial distribution)
 * @param width	Bin width in which counts are averaged.	 */
void LocalNumber::GetHistogram(std::vector<GeometryVector> & result, size_t width) const {
	result.clear();
	if (counts.size() > 0) {
		size_t idx = 0;
		for (idx = 0; idx < counts.size(); idx += width) {
			double x = min + idx + 0.5*(width - 1), x_error = 0.5*width, y_error = 0.0L;
			double y = 0;
			for (int j = 0; j < width; j++) {
				if (idx + j < counts.size()) {
					y += counts[idx + j];
				}
			}
			y /= totalCount;
			y_error = sqrt(y*(1.0 - y) / (this->totalCount - 1.0));
			result.emplace_back(x, y, x_error, y_error);
		}
	}
	else {
		std::cerr << "LocalNumber::GetHistogram()\tThis object is empty.\n";
		return;
	}
}

/** A member function to store data in a binary file. */
void LocalNumber::WriteBinary(std::ostream & ofile) const {
	ofile.write((char *)(&this->min), sizeof(this->min));
	ofile.write((char *)(&this->max), sizeof(this->max));
	ofile.write((char *)(&this->totalCount), sizeof(this->totalCount));
	for (size_t i = 0; i < counts.size(); i++) {
		ofile.write((char *)(&this->counts[i]), sizeof(this->counts[i]));
	}
}
/** A constructor by loading data written by WriteBinary(). */
LocalNumber::LocalNumber(std::istream & ifile) {
	ifile.read((char *)(&this->min), sizeof(this->min));
	ifile.read((char *)(&this->max), sizeof(this->max));
	ifile.read((char *)(&this->totalCount), sizeof(this->totalCount));

	size_t size = this->max - this->min + 1;
	this->counts.resize(size, 0);
	for (size_t i = 0; i < size; i++) {
		ifile.read((char *)(&this->counts[i]), sizeof(this->counts[i]));
	}
	this->CDF_uptodate = false;
}


// ------------------------------------------------------
//	Member functions in BootstrapSetup class
// ------------------------------------------------------
/** \brief Generate LocalNumber object using the prefix and radius.*/
void BootstrapSetup::Load(LocalNumber & data, double R) const {
	std::string FileName(prefix);
	if (R < 0)
		FileName += ".bin";
	else
		FileName += "__" + std::to_string(R) + ".bin";

	std::fstream ifile(FileName.c_str(), std::fstream::in | std::fstream::binary);
	//File exists
	if (ifile.good())
		data = LocalNumber(ifile);
	//doesn't exist
	else
		data = LocalNumber();

	ifile.close();
}

/** \brief Record data in a LocalNumber object in a binary file, if this->SaveData = true*/
void BootstrapSetup::Save(const LocalNumber & LN, double R) const {
	std::string FileName(prefix);
	if (R < 0)
		FileName += ".bin";
	else
		FileName += "__" + std::to_string(R) + ".bin";

	if (this->SaveData) {
		std::fstream ofile(FileName.c_str(), std::fstream::out | std::fstream::binary);
		LN.WriteBinary(ofile);
		ofile.close();
	}
}

/** \brief Compute requested statistics.
 *	Parallelized with OpenMP.*/
void BootstrapSetup::Compute(const LocalNumber & data,
	GeometryVector & Mean, GeometryVector & variance, 
	GeometryVector & skewness, GeometryVector & kurtosis, double R, size_t SampleSize) const {
	std::string FILENAME, output_filename;
	if (R < 0){
		FILENAME = prefix + "__";
		output_filename = output_prefix + "__";
	}
	else{
		FILENAME = prefix + "__" + std::to_string(R) + "__";
		output_filename = output_prefix + "__" + std::to_string(R) + "__";
	}
	data.UpdateCDF();


	if (ComputeVariance || ComputeSkewness || ComputeKurtosis) {
		double sigma2 = data.GetVariance(), mean = data.GetMean();
		double var = 0, skew = 0, kurt = 0;
		double SE_2 = 0, SE_3 = 0, SE_4 = 0;
		double n = 1;
		if((0 < SampleSize) && ((size_t)SampleSize < data.GetTotalCount())){
			n = (double)SampleSize;
		}
		else{
			n = (double)data.GetTotalCount();
		}
		
		double factor_skewness = n*n / ((n - 1.0)*(n - 2.0));
		double factor_kurtosis = n * n*(n*(n - 1.0)*(n - 1.0) + (6.0*n - 9.0)) / ((n - 1.0)*(n - 1.0)*(n - 1.0) *(n*n - 3.0*n + 3.0)),
			const_kurtosis = -3.0 - n * n*n*(6.0*n - 9.0) / ((n - 1.0)*(n - 1.0)*(n - 1.0) *(n*n - 3.0*n + 3.0));

		mean -= data.GetMin();
		//3rd moment
		auto f3 = [&mean](double x)->double {
			double val = x - mean;
			return val * val*val;
		};
		//4th moment
		auto f4 = [&mean](double x)->double {
			double val = x - mean;
			val *= val;	// (x-mean)^2
			val *= val;	// (x-mean)^4
			return val; };

#pragma omp parallel for schedule(guided)
		for (int i = 0; i < BootstrapSize; i++) {
			std::unique_ptr<LocalNumber> ptr = data.GetBootstrapSample(i,SampleSize);
			if (ComputeVariance) {
				//Bootstrap sample variance
				double second = ptr->GetVariance();
#pragma omp atomic
				var += second;
#pragma omp atomic
				SE_2 += second * second;
			}
			if (ComputeSkewness) {
				//Bootstrap sample skewness
				double third = factor_skewness* ptr->ComputeAverageOf(f3) / (sigma2*sqrt(sigma2));
#pragma omp atomic
				skew += third;
#pragma omp atomic
				SE_3 += third * third;
			}
			if (ComputeKurtosis) {
				//Bootstrap sample kurtosis
				double fourth = factor_kurtosis*ptr->ComputeAverageOf(f4) / (sigma2*sigma2) +const_kurtosis;
#pragma omp atomic
				kurt += fourth;
#pragma omp atomic
				SE_4 += fourth * fourth;
			}
		}

		double val, error, sigma = -1;
		if (ComputeMean) {
			val = mean + data.GetMin();
			error = data.GetSE_Mean(val);

			Mean = GeometryVector(val, error);
		}

		if (ComputeVariance) {
			val = var;	error = SE_2;

			val /= BootstrapSize;
			error = (error / BootstrapSize - val * val)*(double)BootstrapSize / (BootstrapSize - 1);
			error = (error > 0) ? sqrt(error) : 0.0;
			variance = GeometryVector(val, error);
		}

		if (ComputeSkewness) {
			val = skew;	error = SE_3;

			val /= BootstrapSize;
			error = (error / BootstrapSize - val * val)*(double)BootstrapSize / (BootstrapSize - 1);
			error = (error > 0) ? sqrt(error) : 0.0;

			skewness = GeometryVector(val, error);
		}

		if (ComputeKurtosis) {
			val = kurt;	error = SE_4;

			val /= BootstrapSize;
			error = (error / BootstrapSize - val * val)*(double)BootstrapSize / (BootstrapSize - 1);
			error = (error > 0) ? sqrt(error) : 0.0;

			kurtosis = GeometryVector(val, error);
		}

		/* Obtain the probability distribution function (PDF) 
		 * and the cummulative distrivution function (CDF) from the post-sampling. */
		if(PrintCDF || PrintPDF){
			std::vector<GeometryVector> CDF, PDF;

//			std::cout << "data: "<< data.GetMin()<< ", " << data.GetMax() <<", "<< data.GetTotalCount()<<"\n";

			/*Initialize CDF*/
			if(PrintCDF){
				for(int i=data.GetMin(); i<data.GetMax()+1; i++){
					CDF.emplace_back((double)i,0.0,0.0,0.0);
				}
			}

			/*Initialize PDF*/
			if(PrintPDF){
				for(int i=data.GetMin(); i<data.GetMax()+1; i++){
					PDF.emplace_back((double)i,0.0,0.0,0.0);
				}
			}

			/*Bootstrap sampling*/
#pragma omp parallel for schedule(guided)
			for (int i = 0; i < BootstrapSize; i++) {
				std::vector<GeometryVector> CDF_temp, PDF_temp;
				std::unique_ptr<LocalNumber> ptr = data.GetBootstrapSample(i,SampleSize);
				
/* #pragma omp critical
{
			std::cout << "ptr: "<< ptr->GetMin()<< ", " << ptr->GetMax() <<", "<< ptr->GetTotalCount()<<"\n";
} */

				if(PrintCDF){
					ptr->UpdateCDF();
					ptr->GetCDF(CDF_temp);

					for (int j=0; j<CDF.size(); j++){
#pragma omp atomic
						CDF[j].x[2] += CDF_temp[j].x[1];
#pragma omp atomic
						CDF[j].x[3] += CDF_temp[j].x[1]*CDF_temp[j].x[1];
						/*save the PDF of the first boostrap sample*/
						if(i==0){
#pragma omp atomic
							CDF[j].x[1] += CDF_temp[j].x[1];
						}
					}
				}

				if(PrintPDF){
					ptr->GetHistogram(PDF_temp, 1);

					for (int j=0; j<PDF.size(); j++){
#pragma omp atomic
						PDF[j].x[2] += PDF_temp[j].x[1];
#pragma omp atomic
						PDF[j].x[3] += PDF_temp[j].x[1]*PDF_temp[j].x[1];
						/*save the PDF of the first boostrap sample*/
						if(i==0){
#pragma omp atomic
							PDF[j].x[1] += PDF_temp[j].x[1];
						}
					}
				}
			}


			/*summarize*/
			if(PrintCDF){
				for(int i=0; i<CDF.size(); i++){
					double mean = CDF[i].x[2]/(double)BootstrapSize;
					double SE_mean = CDF[i].x[3]/(double)BootstrapSize - mean*mean;
					SE_mean = std::max(0.0, SE_mean); /* it prevents the error from being negative due to round-off errors. */
					SE_mean *= ((double)BootstrapSize)/((double)BootstrapSize - 1.0);
					SE_mean = std::sqrt(SE_mean);
					CDF[i].x[2] = SE_mean;
					CDF[i].x[3] = 0;
				}

				std::vector<std::string> label = { 
					"CDF\n",
					"Total number of Samples = " + std::to_string(data.GetTotalCount()) + "\n",
					"Number of bootstrap samples = " + std::to_string(BootstrapSize) + "\n",
					"Size of a bootstrap sample = " + std::to_string(SampleSize) + "\n",
					"N\tCDF(N)\t error in CDF"
				};

				WriteFunction(CDF, (output_filename + "CDF").c_str(), label);
			}

			if(PrintPDF){
				/* Unscaled one */
				for(int i=0; i<PDF.size(); i++){
					double mean = PDF[i].x[2]/(double)BootstrapSize;
					double SE_mean = PDF[i].x[3]/(double)BootstrapSize - mean*mean;
					SE_mean = std::max(0.0, SE_mean); /* it prevents the error from being negative due to round-off errors. */
					SE_mean *= ((double)BootstrapSize/((double)BootstrapSize - 1.0));
					SE_mean = std::sqrt(SE_mean);
					PDF[i].x[2] = SE_mean;
					PDF[i].x[3] = 0;
				}

				std::vector<std::string> label = { 
					"PDF\n",
					"Total number of Samples = " + std::to_string(data.GetTotalCount()) + "\n",
					"Number of bootstrap samples = " + std::to_string(BootstrapSize) + "\n",
					"Size of a bootstrap sample = " + std::to_string(SampleSize) + "\n",
					"N\tPDF(N)\t error in PDF"
				};

				WriteFunction(PDF, (output_filename + "PDF").c_str(), label);

				/* Compute the normalized one. */
				if (ComputeVariance) {
					double mean = data.GetMean(), SE_mean = data.GetSE_Mean(mean);
					double var_ = variance.x[0], SE_var_ = variance.x[1];

#pragma omp parallel for
					for (int i = 0; i < PDF.size(); i++) {
						PDF[i].x[0] = (PDF[i].x[0] - mean) / sqrt(var_);
						PDF[i].x[1] *= sqrt(var_);

						PDF[i].x[2] = (PDF[i].x[2] * PDF[i].x[2] + SE_mean * SE_mean
							+ PDF[i].x[0] * PDF[i].x[0] * SE_var_*SE_var_ / (4.0*var_)) / var_;
						PDF[i].x[2] = sqrt(PDF[i].x[2]);

						PDF[i].x[3] = sqrt(SE_var_*SE_var_ / (4.0*var_*var_)*PDF[i].x[1] * PDF[i].x[1]
							+ var_ * PDF[i].x[3] * PDF[i].x[3]);
					}

					std::vector<std::string> label = { 
						"normalized PDF\n",
						"mean = " + std::to_string(mean) + " +- " + std::to_string(SE_mean) + "\n",
						"var_iance = " + std::to_string(var_) + " +- " + std::to_string(SE_var_) + "\n",
						"(x+<x>)/sigma\tsigma p(x)\t error in x\terror in y"
					};
					WriteFunction(PDF, (output_filename + "PDF_normal").c_str(), label);
				}
			}
		}
	}
}

/** \brief Compute distance metrices from the normal (or Gaussian) distribution.
 *	Parallelized with OpenMP.
	*  @param[in] data		A LocalNumber object that stores the empirical distribution to compute statistics.
	*	@param[in] metrics	metric[i].x[0] = metric_i, metric[i].x[1] = SE_metric_i.
	*	@param[in] R		Window radius.	 
	*	@param[in] SampleSize	If SampleSize = 0, the size of a Bootstrap sample is 
							identical to the original size. */
void BootstrapSetup::ComputeMetrics(const LocalNumber & data, std::vector<GeometryVector> & metrics,
	double R, size_t SampleSize) const{
	
	metrics.clear();
	metrics.resize(2, GeometryVector(0.0,0.0,0.0,0.0));

	std::string FILENAME, output_filename;
	if (R < 0){
		FILENAME = prefix + "__";
		output_filename = output_prefix + "__";
	}
	else{
		FILENAME = prefix + "__" + std::to_string(R) + "__";
		output_filename = output_prefix + "__" + std::to_string(R) + "__";
	}
	data.UpdateCDF();

	if(ComputeL2norm || ComputeKL){
		double L2_norm =0.0L, SE_L2 = 0.0L;
		double KL = 0.0L, SE_KL = 0.0L;

		double sigma_data = data.GetVariance();
		double count = 0.0;
		if(sigma_data >1e-20){

	#pragma omp parallel for schedule(guided)
			for (int i = 0; i < BootstrapSize; i++){
				std::vector<GeometryVector> PDF_test, CDF_test;
				std::vector<GeometryVector> PDF_ref, CDF_ref;	/* Gaussian distribution. */

				std::unique_ptr<LocalNumber> ptr = data.GetBootstrapSample(i,SampleSize);
				ptr->UpdateCDF();

				double mean_bootstrap = ptr->GetMean(), sigma_bootstrap = std::sqrt(ptr->GetVariance());

				if(sigma_bootstrap > 1e-20){
#pragma omp atomic
					count++;

					boost::math::normal_distribution<> gauss_ref(mean_bootstrap, sigma_bootstrap);
					//double norm_factor = 1.0 - cdf(gauss_ref, (double)ptr->GetMin()), //cdf(gauss_ref, (double)ptr->GetMax()) - cdf(gauss_ref, (double)ptr->GetMin()),
					//		norm_const = cdf(gauss_ref, (double)ptr->GetMin());
					int N_min = (int)floor(mean_bootstrap - 5.0 * sigma_bootstrap),
						N_max = (int)ceil(mean_bootstrap + 5.0 * sigma_bootstrap);
					N_min = std::max(N_min, 0);
					N_min = std::min(N_min, ptr->GetMin()); 
					//std::max(0, N_min, ptr->GetMin());
					N_max = std::max(ptr->GetMax(), N_max);
					double norm_factor = 0.0;
					for(int x = N_min; x< N_max+1; x++){
						norm_factor += pdf(gauss_ref, x);
					}


					if(ComputeL2norm){
						/* compute the CDFs. */
						ptr->GetCDF(CDF_test);
						CDF_ref.resize((size_t)(N_max - N_min + 1), GeometryVector(0.0,0.0)) ;
						double CDF_left = cdf(gauss_ref, N_min);
						for(int j = 0; j < CDF_ref.size(); j++){
							CDF_ref[j].x[0] = (double)(N_min + j);
							CDF_ref[j].x[1] = pdf(gauss_ref, CDF_ref[j].x[0]) / norm_factor;
						}
						for (int j = 1; j < CDF_ref.size(); j++){
							CDF_ref[j].x[1] += CDF_ref[j-1].x[1];
						}
						double val = Compute_L2_norm(CDF_test, CDF_ref, 1.0 / sigma_bootstrap);

						// if(i==0){
						// 	char name_temp[200];
						// 	sprintf(name_temp, "Test__%0.4f__CDF", R);
						// 	WriteFunction(CDF_ref, name_temp);
						// }

		#pragma omp atomic
						L2_norm += val;
		#pragma omp atomic
						SE_L2 += val * val;
					}

					if(ComputeKL){
						/* compute the PDFs */
						ptr->GetHistogram(PDF_test, 1);
						PDF_ref.resize((size_t)(N_max - N_min + 1), GeometryVector(0.0,0.0,0.0,0.0)) ;
						for(int j = 0; j < PDF_ref.size(); j++){
							PDF_ref[j].x[0] = (double)(N_min + j);
							PDF_ref[j].x[1] = pdf(gauss_ref, PDF_ref[j].x[0]) / norm_factor;
						}
						double val = Compute_KullbackLeibler_distance(PDF_test, PDF_ref);

						// if(i==0){
						// 	char name_temp[200];
						// 	sprintf(name_temp, "Test__%0.4f__PDF", R);
						// 	WriteFunction(PDF_ref, name_temp);
						// }

		#pragma omp atomic
						KL += val;
		#pragma omp atomic
						SE_KL += val * val;
					}

				}

			}
		}
		else{
			SE_L2 = -1;
			SE_KL = -1;
		}



		if(ComputeL2norm){
			if(SE_L2 > 0 && count > 1.0){
				L2_norm /= count;	
				SE_L2 = SE_L2 / count - L2_norm * L2_norm;
				SE_L2 = std::max(0.0, SE_L2);
				SE_L2 = std::sqrt(SE_L2 / (count - 1.0)); 
			}

			metrics[0].x[0] = L2_norm;
			metrics[0].x[1] = SE_L2;
		}

		if(ComputeKL){
			if (SE_KL > 0 && count > 1.0){
				KL /= count;
				SE_KL = SE_KL / count - KL * KL;
				SE_KL = std::max(0.0, SE_KL);
				SE_KL = std::sqrt(SE_KL / (count - 1.0)); 
			}

			metrics[1].x[0] = KL;
			metrics[1].x[1] = SE_KL;
		}

	}
}


/** \brief	A function to give uncorrelated sampling points.*/
std::vector<GeometryVector> GetSamplingPoints(DimensionType dim, size_t NumPoints, int seed) {
	RandomGenerator rng(seed);
	GeometryVector pos; pos.SetDimension(dim);
	std::vector<GeometryVector> result;	result.reserve(NumPoints);
	for (size_t i = 0; i < NumPoints; i++) {
		for (DimensionType d = 0; d < dim; d++) {
			pos.x[d] = rng.RandomDouble();
		}
		result.push_back(pos);
	}
	return result;
}



// ------------------------------------------------------
//		Functions to compute local number variance (Spherical windows).
// ------------------------------------------------------

/** \brief Measure the point numbers within a spherical window of a given maximal radius and a given center.
 *	This function searches all points inside a sphere of a maximum radius,
and compute their distances from the sphere center, which tell the radius of a sphere to which the particle belongs.

 *	@param[in] CurrConfig	A Configuration object.
 *	@parma[out] results		A table of N(R), the number of particles inside spheres.
				results[i] = N(R) for R = dR * (i+1)
 *	@param[in] Rmax			A maximal radius up to which N(R) is computed
 *	@param[in] dR			A bin width.
 *	@param[in] Center_RelativeCoordinates	A center of the window in relative Coordinates. */
inline void GetN(const Configuration & CurrConfig, std::vector<size_t> &results,
	double Rmax, double dR, const GeometryVector & Center_RelativeCoordinates) {

	//Confirm dimensionality.
	if (Center_RelativeCoordinates.Dimension == CurrConfig.GetDimension()) {
		//Initialize results.
		size_t Bin = (size_t)ceil(Rmax / dR);
		results.resize(Bin, 0.0L);

		//Search all points that |x0 - r | <Rmax.
		CurrConfig.IterateThroughNeighbors(Center_RelativeCoordinates, Rmax,
			[&results, &dR, &Bin](const GeometryVector &shift, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void
		{
			//results[i] = # of points between dR*i and dR*(i+1).
			size_t index = (size_t)floor(sqrt(shift.Modulus2()) / dR);
			if (index < Bin)
				results[index]++;
		});
		for (size_t i = 1; i < results.size(); i++) {
			results[i] += results[i - 1];
		}
	}
	else {
		std::cerr << "NumberVariance::GetN  : inconsistent dimensionality. \n";
	}
}





/** \brief Measure the Local number variance associated with spherical windows.
*	This function is specialized to efficiently compute the local number variance.
*	@param[in] GetConfig	A lambda function to generate Configuration objects.
*	@param[in] NumConfigs	The number of Configurations.
*	@param[out] Result		A list of GeometryVector objects that describe the local number variances.
				Result[i].x[0] = R
				Result[i].x[1] = sigma_N ^2(R)
				Result[i].x[2] = SE in sigma_N ^2(R), which is computed under the assumption that N(R) follows the central limit theorem.
*	@param[in] Rmax			The largest window radius
*	@param[in] dR			R-resolution
*	@param[in] WindowCenters	GeometryVector objects that describe the centers of windows in Relative coordinates.
								Users can prescribe the centers of sampling windows through this argumement,
								but also use default centers that are binomial point pattern of "NumWindowCenters" particles.  */
void MC_nv(const std::function<Configuration(size_t i)> & GetConfig, size_t NumConfigs,
	std::vector<GeometryVector> & Result, double Rmax, double dR, const std::vector<GeometryVector> & WindowCenters) {

	Configuration CurrConfig = GetConfig(0);
	DimensionType dim = CurrConfig.GetDimension();

	//Centers of sampling windows.
	std::vector<GeometryVector> x0(WindowCenters);

	//Generate Bins
	//Bins[i].x[0] = sum of N(R)
	//Bins[i].x[1] = sum of N(R)^2
	//Bins[i].x[2] = sum of variance
	//Bins[i].x[2] = sum of variance^2

	size_t NumBins = (size_t)ceil(Rmax / dR);
	std::vector<GeometryVector> Bins(NumBins, GeometryVector(0.0L, 0.0L, 0.0L, 0.0L));

	std::cout << "Computing the local number variance \n";
	progress_display pd(NumConfigs);
	//For each configuration,
	for (size_t i = 0; i < NumConfigs; i++) {
		if (i > 0)
			CurrConfig = GetConfig(i);
		if (CurrConfig.NumParticle() == 0) {
			std::cerr << "The particle number should be positive!\n";
			return;
		}
		//Bins for the configuration i;
		std::vector<GeometryVector> Bins_i(NumBins, GeometryVector(0.0L, 0.0L, 0.0L, 0.0L));

		//Adjust cell sizes of CurrConfig
		//to make IterateThroughNeighbors() thread-safe.
		CurrConfig.PrepareIterateThroughNeighbors(Rmax);

#pragma omp parallel
		{
			//Copy "CurrConfig" in each thread to ensure that there is no conflict
			//in using IterateThroughNeighbors(...) function.
			Configuration c(CurrConfig);

#pragma omp for schedule(guided)
			for (int num_samp = 0; num_samp < x0.size(); num_samp++) {
				std::vector<size_t> NR;
				GetN(c, NR, Rmax, dR, x0[num_samp]);

				for (size_t j = 0; j < Bins_i.size(); j++) {
#pragma omp atomic
					Bins_i[j].x[0] += NR[j];
#pragma omp atomic
					Bins_i[j].x[1] += NR[j] * NR[j];
				}
			}
		}

		//Compute variance for this configuraiton
		for (size_t idx_bin = 0; idx_bin < Bins_i.size(); idx_bin++) {
			Bins_i[idx_bin].x[0] /= (x0.size());
			Bins_i[idx_bin].x[1] /= (x0.size());
			double variance = (Bins_i[idx_bin].x[1] - Bins_i[idx_bin].x[0] * Bins_i[idx_bin].x[0])*(x0.size() / (x0.size() - 1.0));
			variance = (variance > 0) ? variance : 0.0;	//check
			Bins_i[idx_bin].x[2] = variance;
			Bins_i[idx_bin].x[3] = variance * variance;

			//Sum Bins_i to Bins
			Bins[idx_bin] = Bins[idx_bin] + 1.0 / (i + 1) *(Bins_i[idx_bin] - Bins[idx_bin]);
		}
		//Increase a count in the display
		pd++;
	}

	//Summarize data
	double NumTotalSampling = NumConfigs * WindowCenters.size();
	Result.resize(Bins.size(), GeometryVector(0, 0, 0, 0));
	for (size_t i = 0; i < Bins.size(); i++) {
		double R = dR * (i + 1);
		Result[i].x[0] = R;

		// Compute the mean and variance
		double mean = Bins[i].x[0];
		double variance = (Bins[i].x[1] - mean * mean)*(NumTotalSampling / (NumTotalSampling - 1.0));
		Result[i].x[1] = (variance > 0) ? variance : 0.0;

		// Compute error
		double mean_var = Bins[i].x[2];
		//Standard error in sample variances.
		double error1 = Bins[i].x[3] - mean_var * mean_var;
		error1 = (error1 > 0) ? sqrt(error1*NumTotalSampling / (NumTotalSampling - 1.0)) : 0.0;
		//Estimated standard error if N(R) follows Gaussian distribution
		double error2 = sqrt(2.0 / (NumTotalSampling - 1.0) * variance);
		Result[i].x[2] = std::max(error1, error2);	//Choose the larger one.

	}

}



/** \brief Compute Local number variance and other statistics by taking advantages of Bootstrap method.
*
*	@param[in] GetConfig	A lambda function to construct Configuration objects.
*	@param[in] NumConfigs	The number of configurations.
*	@param[out] Results		A table of local number variance. @see MC_Nv()
*	@param[in] Rmax			Largest window radius.
*	@param[in] dR			Bin width of window radius
*	@param[in] WHAT2DO		A BootstrapSetup class that stores information about what to compute and aids certain bootstrap computations.
*	@param[in] NumWindowCenters	The number of sampling windows.
*	@param[in] see			A seed for a random number generator.
*	@param[in] WindowsCenters	A list of sampling centers in relative coordinates.
Users can prescribe the centers of sampling windows through this argumement,
but also use default centers that are binomial point pattern of "NumWindowCenters" particles.  */
void MC_nv_distribution(const std::function<Configuration(size_t i)> & GetConfig, size_t NumConfigs,
	std::vector<GeometryVector> & Result, double Rmax, double dR, const BootstrapSetup & WHAT2DO,
	const std::vector<GeometryVector> & WindowCenters) {

	Configuration CurrConfig = GetConfig(0);
	DimensionType dim = CurrConfig.GetDimension();

	//Centers of sampling windows.
	std::vector<GeometryVector> x0(WindowCenters);

	//Generate LocalNumber objects.
	//Bins[i] will store distribution of N(dR*(i+1))
	//Load previous data if they exist.
	std::vector<LocalNumber> Bins;
	{
		std::cout << "load data: ";
		size_t NumBins = (size_t)floor(Rmax / dR);
		Bins.resize(NumBins);
		for (size_t i = 0; i < NumBins; i++) {
			double R = dR * (i + 1);
			WHAT2DO.Load(Bins[i], R);

			if (Bins[i].GetTotalCount() > 0)
				std::cout << "."; 
		}
		std::cout << "\n";
	}

	std::cout << "Computing the local number variance\n";
	progress_display pd(NumConfigs);

	//For each configuration,
	for (size_t i = 0; i < NumConfigs; i++) {
		if (i > 0)
			CurrConfig = GetConfig(i);

		//Adjust cell sizes of CurrConfig
		//to make IterateThroughNeighbors() thread-safe.
		CurrConfig.PrepareIterateThroughNeighbors(Rmax);

#pragma omp parallel for schedule( guided)
		for (int idx_x0 = 0; idx_x0 < x0.size(); idx_x0++) {
			std::vector<size_t> NR;
			GetN(CurrConfig, NR, Rmax, dR, x0[idx_x0]);

			for (size_t j = 0; j < Bins.size(); j++) {
#pragma omp critical
				{
					Bins[j].report(NR[j]);
				}
			}

		}
#pragma omp critical
		{
			pd++;
		}
	}

	//Summarize data
	Result.resize(Bins.size(), GeometryVector(0, 0, 0, 0));
	std::vector<GeometryVector> mean, skewness, kurtosis;
	for (size_t i = 0; i < Bins.size(); i++) {
		double R = dR * (i + 1);
		Result[i].x[0] = R;

		//First save data with a prescribed directory.
		WHAT2DO.Save(Bins[i], R);
		{
			GeometryVector av(0), var(0), skew(0), kurt(0);
			//Do computations.
			//The following function is parallized.
			WHAT2DO.Compute(Bins[i], av, var, skew, kurt, R);
			//If var is computed.
			if (var.Dimension == 2) {
				Result[i].x[1] = var.x[0];
				Result[i].x[2] = var.x[1];
			}
			//If not, we simply compute it...
			else {
				Result[i].x[1] = Bins[i].GetVariance();
			}

			//If mean is computed
			if (av.Dimension == 2) {
				mean.emplace_back(R, av.x[0], av.x[1]);
			}

			//If skew is computed
			if (skew.Dimension == 2) {
				skewness.emplace_back(R, skew.x[0], skew.x[1]);
			}

			//If kurtosis is computed
			if (kurt.Dimension == 2) {
				kurtosis.emplace_back(R, kurt.x[0], kurt.x[1]);
			}
		}
	}
	if (WHAT2DO.ComputeMean) {
		std::vector<std::string> label = {
			"unbaised estimaters of mean\n",
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize) + "\n",
			"R\tmean\tSE_skewness"
		};
		WriteFunction(mean, (WHAT2DO.prefix + "__mean").c_str(), label);
	}
	if (WHAT2DO.ComputeSkewness) {
		std::vector<std::string> label = { 
			"unbaised estimaters of skewness\n", 
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize)+"\n",
			"R\tskewness\tSE_skewness"
		};
		WriteFunction(skewness, (WHAT2DO.prefix + "__skewness").c_str(), label);
	}
	if (WHAT2DO.ComputeKurtosis) {
		std::vector<std::string> label = {
			"unbiased estimator of kurtosis\n",
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize) + "\n",
			"R\texcess kurtosis\tSE_kurtosis"
		};
		WriteFunction(kurtosis, (WHAT2DO.prefix + "__kurtosis").c_str(),label);
	}
	
}


/** \brief Compute Local number variance and other statistics by taking advantages of Bootstrap method.
 *	This function is parallelized in the for loops of particles.
 *	This function implements a post sampling method to reduce the finite-size effects
 *	in the void MC_nv_distribution(...) function.
 *	@param[in] GetConfig	A lambda function to construct Configuration objects.
 *	@param[in] NumConfigs	The number of configurations.
 *	@param[out] Results		A table of local number variance. @see MC_Nv()
 *	@param[in] Rmax			Largest window radius.
 *	@param[in] dR			Bin width of window radius
 *	@param[in] WHAT2DO		A BootstrapSetup class that stores information about what to compute and aids certain bootstrap computations.
 *	@param[in] WindowsCenters	A list of sampling centers in relative coordinates.
								Users can prescribe the centers of sampling windows through this argumement.  */
void MC_nv_distribution_PostSampling(const std::function<Configuration(size_t i)> & GetConfig, size_t NumConfigs,
	std::vector<GeometryVector> & Result, double Rmax, double dR, const BootstrapSetup & WHAT2DO, const std::vector<GeometryVector> & WindowCenters){

	Configuration CurrConfig = GetConfig(0);
	int dim = CurrConfig.GetDimension();

	double tot_num_prts = 0;	/*total number of particles.*/ 
	double rho_overall = 0;		/*overall number density. */
	//Centers of sampling windows.
	std::vector<GeometryVector> x0(WindowCenters);

	//Generate LocalNumber objects.
	//Bins[i] will store distribution of N(dR*(i+1))
	//Load previous data if they exist.
	std::vector<LocalNumber> Bins;
	{
		std::cout << "load data: ";
		size_t NumBins = (size_t)floor(Rmax / dR);
		Bins.resize(NumBins);
		for (size_t i = 0; i < NumBins; i++) {
			double R = dR * (i + 1);
			WHAT2DO.Load(Bins[i], R);

			if (Bins[i].GetTotalCount() > 0)
				std::cout << "."; 
		}
		std::cout << "\n";
	}

	std::cout << "Computing the local number variance\n";
	progress_display pd(NumConfigs);

	//For each configuration,
	for (size_t i = 0; i < NumConfigs; i++) {
		if (i > 0)
			CurrConfig = GetConfig(i);

		//Adjust cell sizes of CurrConfig
		//to make IterateThroughNeighbors() thread-safe.
		CurrConfig.PrepareIterateThroughNeighbors(Rmax);

		tot_num_prts += (double)CurrConfig.NumParticle();
		rho_overall += (double)CurrConfig.NumParticle()/CurrConfig.PeriodicVolume();
#pragma omp parallel for schedule( guided)
		for (int idx_x0 = 0; idx_x0 < x0.size(); idx_x0++) {
			std::vector<size_t> NR;
			GetN(CurrConfig, NR, Rmax, dR, x0[idx_x0]);

			for (size_t j = 0; j < Bins.size(); j++) {
#pragma omp critical
				{
					Bins[j].report(NR[j]);
				}
			}

		}
#pragma omp critical
		{
			pd++;
		}
	}

	tot_num_prts /= (double)NumConfigs;
	rho_overall /= (double)NumConfigs;

	//Post-process data
	Result.resize(Bins.size(), GeometryVector(0, 0, 0, 0));
	std::vector<GeometryVector> mean, skewness, kurtosis;
	std::vector<GeometryVector> L2, KL;
	for (size_t i = 0; i < Bins.size(); i++) {
		double R = dR * (i + 1);
		Result[i].x[0] = R;

		/*------------------------------------------
			Determine a proper number of windows
			N_window <= tot_num_prt * NumConfig / (2*rho_overall* v1(R))
		*------------------------------------------*/
		size_t N_window = std::max(
			(size_t)round(tot_num_prts * NumConfigs /(2.0*rho_overall * HyperSphere_Volume(dim,R))), 
			(size_t)1);
		if(N_window > Bins[i].GetTotalCount())
			N_window = 0;
		//First save data in a prescribed directory.
		WHAT2DO.Save(Bins[i], R);

		int N_window_per_config;
		/* ---- the following line disables the post-sampling !! ---- */
		//N_window = 0; 
		/* ---- turned on the bootstrap sampling ---- */
		N_window = Bins[i].GetTotalCount();

		if (N_window == 0)
			N_window_per_config = (int)Bins[i].GetTotalCount()/NumConfigs;
		else 
			N_window_per_config = N_window / NumConfigs;


		/* Compute moment-related statistics (variance, skewness, ex. kurtosis).*/
		{
			GeometryVector av(0), var(0), skew(0), kurt(0);
			//Do computations.
			//The following function is parallized.
			if(N_window == 0)
				WHAT2DO.Compute(Bins[i], av, var, skew, kurt, R);
			else
			{
				WHAT2DO.Compute(Bins[i], av, var, skew, kurt, R, N_window);
			}
			
			
			//If var is computed.
			if (var.Dimension == 2) {
				Result[i].x[1] = var.x[0];
				Result[i].x[2] = var.x[1];
				Result[i].x[3] = N_window_per_config;
			}
			//If not, we simply compute it...
			else {
				Result[i].x[1] = Bins[i].GetVariance();
			}

			//If mean is computed
			if (av.Dimension == 2) {
				mean.emplace_back(R, av.x[0], av.x[1]);
			}

			//If skew is computed
			if (skew.Dimension == 2) {
				skewness.emplace_back(R, skew.x[0], skew.x[1], N_window_per_config);
			}

			//If kurtosis is computed
			if (kurt.Dimension == 2) {
				kurtosis.emplace_back(R, kurt.x[0], kurt.x[1], N_window_per_config);
			}
		}
		/* distance metrics are computed. */
		{	
			std::vector<GeometryVector> metrics;
			//Do computations.
			//The following function is parallized.
			if(N_window == 0)
				WHAT2DO.ComputeMetrics(Bins[i], metrics, R);
			else
				WHAT2DO.ComputeMetrics(Bins[i], metrics, R, N_window);
			
			//If L2-norm is computed.
			if (WHAT2DO.ComputeL2norm) {
				if(metrics[0].x[1] > 0){
					L2.emplace_back(R, metrics[0].x[0], metrics[0].x[1], N_window_per_config);
				}
			}
			//If Kullback-Leibler distance is computed. 
			if (WHAT2DO.ComputeKL){
				if(metrics[1].x[1] > 0){
					KL.emplace_back(R, metrics[1].x[0], metrics[1].x[1], N_window_per_config);
				}
			}
		}
	}
	if (WHAT2DO.ComputeMean) {
		std::vector<std::string> label = {
			"unbaised estimaters of mean\n",
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize) + "\n",
			"R\tmean\tSE_mean"
		};
		WriteFunction(mean, (WHAT2DO.output_prefix + "__mean").c_str(), label);
	}
	if (WHAT2DO.ComputeVariance) {
		std::vector<std::string> label = { 
			"Post-sampling with a fixed size\n",
			"unbaised estimaters of local number variance\n", 
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize)+"\n",
			"R\tsigma_N^2(R)\tSE_sigma_N^2(R)\tN_window/config"
		};
		WriteFunction(Result, (WHAT2DO.output_prefix + "__LNV2").c_str(), label);
	}
	if (WHAT2DO.ComputeSkewness) {
		std::vector<std::string> label = { 
			"Post-sampling with a fixed size\n",
			"unbaised estimaters of skewness\n", 
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize)+"\n",
			"R\tskewness\tSE_skewness\tN_window/config"
		};
		WriteFunction(skewness, (WHAT2DO.output_prefix + "__skewness").c_str(), label);
	}
	if (WHAT2DO.ComputeKurtosis) {
		std::vector<std::string> label = {
			"Post-sampling with a fixed size\n",
			"unbiased estimator of excess kurtosis\n",
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize) + "\n",
			"R\texcess kurtosis\tSE_kurtosis\tN_window/config"
		};
		WriteFunction(kurtosis, (WHAT2DO.output_prefix + "__kurtosis").c_str(),label);
	}
	if (WHAT2DO.ComputeL2norm) {
		std::vector<std::string> label = {
			"Post-sampling with a fixed size\n",
			"L2-norm of the CDF of bootstrap samples w.r.t. Gaussian distributions.\n",
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize) + "\n",
			"Standard Error (SE) is estimated from bootstrap samples. \n",
			"R\tL2-norm\tSE_L2-norm\tN_window/config"
		};
		WriteFunction(L2, (WHAT2DO.output_prefix + "__L2").c_str(),label);
	}	
	if (WHAT2DO.ComputeKL) {
		std::vector<std::string> label = {
			"Post-sampling with a fixed size\n",
			"Kullback-Leibler distance (relative entropy) of PDF w.r.t. Gaussian distribution.\n",
			"total samples = " + std::to_string(Bins[0].GetTotalCount()) + "\n",
			"number of bootstrap samples = " + std::to_string(WHAT2DO.BootstrapSize) + "\n",
			"Standard Error (SE) is estimated from bootstrap samples. \n",
			"R\tK.L. \tSE_K.L.\tN_window/config"
		};
		WriteFunction(KL, (WHAT2DO.output_prefix + "__KL").c_str(),label);
	}

}

double Compute_L2_norm(const std::vector<GeometryVector> & CDF_test, 
					const std::vector<GeometryVector> & CDF_ref, double dx) {
	if(CDF_test.size() <= CDF_ref.size()){
		double sum = 0.0L, delta = 0.0;
		
		int j=0;
		for(int i=0; i< CDF_ref.size(); i++){
			
			if( (j < CDF_test.size()) && 
				(abs(CDF_test[j].x[0] - CDF_ref[i].x[0]) < 1e-3) ){
				delta = CDF_test[j].x[1] - CDF_ref[i].x[1];
				sum += delta * delta;
				j++;
			}
		}
		sum = sum*dx;
		return sum;
	}
	else{
		std::cerr<<"Compute_L2_norm::The size of reference distribution must be larger than those of the test distributions.\n";
		return 0;
	}
}

double Compute_KullbackLeibler_distance(const std::vector<GeometryVector> & PDF_test, 
	const std::vector<GeometryVector> & PDF_ref, double dx){

	if(PDF_test.size() <= PDF_ref.size()){
		double KL = 0.0L, temp = 0.0L;
		
		int j=0;
		for(int i=0; i < PDF_ref.size(); i++){

			if( (j < PDF_test.size()) && 
				(abs(PDF_test[j].x[0] - PDF_ref[i].x[0]) < 1e-3)  ){
				
				if(PDF_test[j].x[1] > 1e-10){
					temp = PDF_test[j].x[1] * log(PDF_test[j].x[1] / PDF_ref[i].x[1]);
					KL += temp;
				}
				j++;
			}		
		}
		return KL;
	}
	else{
		std::cerr<<"Compute_KullbackLeibler_distance::The size of reference distribution must be larger than those of the test distributions.\n";
		return 0;
	}	

}
