/**	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	:	April. 2019 */

/** \file NumberVariance.h 
 *	\brief Header file to compute the local number variance from Configuration objects.
*/


#ifndef __NUMBERVARIANCE_H__
#define __NUMBERVARIANCE_H__


#include <vector>
#include <omp.h>
#include <cassert>
#include <functional>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>

/* header files in $(cores) */
#include <GeometryVector.h>
#include <PeriodicCellList.h>
#include <RandomGenerator.h>


/// @brief Function to compute the smallest radius for a spherical 'window' to contain a displacement vector from the origin. This function allows one to implement arbitrarily shaped windows.
/// @param delta_x displacement vector from the window center to the position of a particle.
/// @return the defined distance.
inline double radius_in_Sphere (const GeometryVector & delta_x){
	return std::sqrt(delta_x.Modulus2());
}


/// @brief Compute the superdisk window radius in which a point (relative to the window center) is located. i.e., x^(2*p) + y^(2*p) + z^(2*p) = b^(2*p) 
/// @param delta_x 
/// @param p
/// @return b
inline double radius_in_Superdisk (const GeometryVector & delta_x, double p){
	double distance2 = 0.0;
	for (int i=0; i<delta_x.Dimension; i++){
		distance2 += pow(delta_x.x[i], 2.*p);
	}
	return exp(log(distance2) / (2.*p));
}




// ------------------------------------------------------
//		Functions to define sampling points.
// ------------------------------------------------------

/** \brief	A function to give uncorrelated sampling points.
	@param dim			Space dimensions
	@param NumPoints	The number of points
	@param centers  	Window centers
	@param rng			A random generator.
	*/
void GetSamplingPoints(DimensionType dim, size_t NumPoints, std::vector<GeometryVector> & centers , RandomGenerator & rng);


/** \brief Measure the point numbers within a given window of a given maximal radius and a given center.
 *	This function searches all points inside a sphere of a maximum radius, 
	and computes their distances from the sphere center, which tell the radius of a sphere to which the particle belongs.

 *	@param[in] config	A Configuration object.
 *	@param[out] results		A table of N(R), the number of particles inside spheres.
							results[i] = N(R) for R = dR * (i+1)
 *	@param[in] Rmax			A maximal radius up to which N(R) is computed
 *	@param[in] dR			A bin width.
 *	@param[in] x0_RelativeCoordinates	A center of the window in relative Coordinates. 
 *	@param[in] WindowSizeFunct	function from a vector to float, which define the size of window. By default, it is Eucliean distance, leading to spherical window.
 */
void GetN(const Configuration & config, std::vector<double> &results,
	double Rmax, double dR, const GeometryVector & x0_RelativeCoordinates, const std::function<double(const GeometryVector &)> WindowSizeFunct = radius_in_Sphere);


/// @brief Monte Carlo calculation of the local number variance of a single configuration. Window centers are given.
/// @param config 	
/// @param results 	results[i] = {R, variance}
/// @param Rmax 	Largest window size
/// @param dR 		Increment in window size
/// @param window_centers 	Specified window centers.
/// @param WindowSizeFunct 
void MC_NV_single(const Configuration &  config, std::vector<GeometryVector> & results, double Rmax, double dR, const std::vector<GeometryVector> &  window_centers, const std::function<double(const GeometryVector & )> WindowSizeFunct = radius_in_Sphere);

void MC_NV_single(const Configuration & config, std::vector<GeometryVector> &results, double Rmax, double dR, size_t num_window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct = radius_in_Sphere, int seed = 0);


/** \brief Measure the Local number variance associated with spherical windows. 
 *	This function is specialized to efficiently compute the local number variance.
 *	@param[in] GetConfig	A lambda function to generate Configuration objects.
 *	@param[in] NumConfigs	The number of Configurations.
 *	@param[out] Result		A list of GeometryVector objects that describe the local number variances.
						Result[i].x = {R, variance, SE in variance}
						The SE is computed under the assumption that N(R) follows the central limit theorem.
 *	@param[in] Rmax			The largest window radius
 *	@param[in] dR			Resolution in R
 *	@param[in] WindowCenters	GeometryVector objects that describe the centers of windows in Relative coordinates.
					Users can prescribe the centers of sampling windows through this argumement,
					but also use default centers that are binomial point pattern of "NumWindowCenters" particles.  
*	@param[in] WindowSizeFunct 				
					*/
void MC_NV_Ensemble(const std::function<Configuration(size_t i)> & GetConfigsFunction, size_t NumConfigs, 
	std::vector<GeometryVector> & Result, double Rmax, double dR, const std::vector<GeometryVector> & WindowCenters, const std::function<double(const GeometryVector &)> WindowSizeFunct = radius_in_Sphere);

/** \brief Measure the Local number variance associated with the prescribed windows. 
 *	This function is specialized to efficiently compute the local number variance.
 *	@param[in] GetConfig	A lambda function to generate Configuration objects.
 *	@param[in] NumConfigs	The number of Configurations.
 *	@param[out] Result		A list of GeometryVector objects that describe the local number variances.
						Result[i].x = {R, variance, SE in variance}
						The SE is computed under the assumption that N(R) follows the central limit theorem.
 *	@param[in] Rmax			The largest window radius
 *	@param[in] dR			Resolution in R
 *	@param[in] WindowCenters	GeometryVector objects that describe the centers of windows in Relative coordinates.
					Users can prescribe the centers of sampling windows through this argumement,
					but also use default centers that are binomial point pattern of "NumWindowCenters" particles.  
*	@param[in] WindowSizeFunct 				
*	@param[in] change_centers	Use different window centers in each configuration. 				
					*/
void MC_NV_Ensemble(const std::function<Configuration(size_t i)> & GetConfigsFunction, size_t NumConfigs, 
	std::vector<GeometryVector> & Result, double Rmax, double dR, size_t num_window_centers, const std::function<double(const GeometryVector &)> WindowSizeFunct = radius_in_Sphere, int seed=0, bool change_centers = false);


//TODO: Add functions for random orientations.
/** \brief Measure the point numbers within a given window of a given maximal radius, a given center, and a random orientation
 *	This function searches all points inside a sphere of a maximum radius, 
	and computes their distances from the sphere center, which tell the radius of a sphere to which the particle belongs.

 *	@param[in] config	A Configuration object.
 *	@param[out] results		A table of N(R), the number of particles inside spheres.
							results[i] = N(R) for R = dR * (i+1)
 *	@param[in] Rmax			A maximal radius up to which N(R) is computed
 *	@param[in] dR			A bin width.
 *	@param[in] x0_RelativeCoordinates	A center of the window in relative Coordinates. 
 *	@param[in] WindowSizeFunct	function from a vector to float, which define the size of window. By default, it is Eucliean distance, leading to spherical window.
//  */
// void GetN_orientation(const Configuration & config, std::vector<double> &results,
// 	double Rmax, double dR, const GeometryVector & x0_RelativeCoordinates, const std::function<double(const GeometryVector &)> WindowSizeFunct = radius_in_Sphere, RandomGenerator & rng);


/** \brief Measure the Local number variance associated with the prescribed windows and their random orientations.
 *	This function is specialized to efficiently compute the local number variance.
 *	@param[in] GetConfig	A lambda function to generate Configuration objects.
 *	@param[in] NumConfigs	The number of Configurations.
 *	@param[out] Result		A list of GeometryVector objects that describe the local number variances.
						Result[i].x = {R, variance, SE in variance}
						The SE is computed under the assumption that N(R) follows the central limit theorem.
 *	@param[in] Rmax			The largest window radius
 *	@param[in] dR			Resolution in R
 *	@param[in] WindowCenters	GeometryVector objects that describe the centers of windows in Relative coordinates.
					Users can prescribe the centers of sampling windows through this argumement,
					but also use default centers that are binomial point pattern of "NumWindowCenters" particles.  
*	@param[in] WindowSizeFunct 				
*	@param[in] change_centers	Use different window centers in each configuration. 				
					*/
void MC_NV_Ensemble_Orientational(const std::function<Configuration(size_t i)> & GetConfigsFunction, size_t NumConfigs, 
	std::vector<GeometryVector> & Result, double Rmax, double dR, const std::vector<GeometryVector> & WindowCenters, const std::function<double(const GeometryVector &)> WindowSizeFunct = radius_in_Sphere);


#endif