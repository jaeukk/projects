/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: April 2021 */


/** \file FirstPassageTime.h
  * \brief Header file for a FirstPassageTime simulation
  * for two-phase packings.
  * */

#ifndef FIRSTPASSAGETIME_H
#define FIRSTPASSAGETIME_H

#include "../NonsphericalParticles/Dispersions.h"
#include "../PeriodicCellList.h"
#include "../GeometryVector.h"
#include "../RandomGenerator.h"
#include "../etc.h"

#include <omp.h>

//TODO: Superconducting phase
//TODO: Insulting phase

/** \struct A class to store simple data for trajectories. */
struct Trajectory{
	DimensionType dimension = 1;
	double sigma_host = 0.0, sigma_particle = 0.0;
	double R2_host = 0.0, R2_particle = 0.0;
	double R2_s = 0.0, tau_s = 0.0;	/* data near the interface*/
	double tau_total = 0.0;

	GeometryVector pos_curr;	/* test particle position in the relative coordinates. */
	GeometryVector delta_pos;   /* a displacement from the initial position in the Cartesian */

	Trajectory(){
		sigma_host = 0.0; sigma_particle = 0.0;
		R2_host = 0.0;	R2_particle = 0.0;
		R2_s = 0.0; tau_s = 0.0;
	}

	Trajectory(DimensionType d, double sigma1, double sigma2){
		dimension = d;
		sigma_host = sigma1;
		sigma_particle = sigma2;

		pos_curr = GeometryVector(static_cast<int> (d));
		delta_pos = GeometryVector(static_cast<int> (d));
	}

	Trajectory(const Trajectory & src);

	void ReadLog(std::istream & ifile);

	inline double GetX2() const{
		return delta_pos.Modulus2();
	}

	inline void Update_R2_host(double R2_h, double tau_assigned = 0.0){
		
		R2_host += R2_h;
		if (sigma_host == 0.0){
			tau_total += tau_assigned;
		}
		else{
			tau_total += R2_h / (2.0*this->dimension * sigma_host);
		}
	}

	inline void Update_R2_particle(double R2_p, double tau_assigned = 0.0){
		R2_particle += R2_p;
		if (sigma_particle == 0.0){
			tau_total += tau_assigned;
		}
		else{
			tau_total += R2_p / (2.0*this->dimension * sigma_particle);
		}
	}

	inline void Update_R2_s(double R2_surface, double tau_surface){
		this->R2_s += R2_surface;
		this->tau_s += tau_surface;
		this->tau_total += tau_surface;
	}

	inline double GetTimeTotal() const{
		return this->tau_total;
	}
};

std::ostream & operator << (std::ostream & os, const Trajectory & t);

struct FirstPassageTimeParameters{
	double sigma_host, sigma_particle;
	double alpha; 
	
	size_t dimension;
	size_t num_test_prts;
	size_t num_threads;

	FirstPassageTimeParameters(double sigma_host, double sigma_particle, size_t Dimension){
		this->sigma_host = sigma_host; this->sigma_particle = sigma_particle;
		this->dimension = Dimension;
		
		if(sigma_host*sigma_particle > 0.0){
			this->alpha = sigma_particle / sigma_host;}
		else{
			this->alpha = 0.0;
		}
	}

	bool ContinueTrajectory(Trajectory & traj);

	void WriteLogs(const std::string & name, const std::vector<Trajectory> & trajs);
	bool ReadLogs(const std::string & name, std::vector<Trajectory> & trajs);

/** @brief The probability that a test particle goes to the host phase. 
 * The probability that a test particle in phase 1 (\theta = \pi/2) or phase 2 (\theta = 3\pi/2) hits \partial \Omega_1 in 2D for the first time.
 * It corresponds to p_1(r, \theta = \pi/2) or p_1(r, \theta = 3\pi/2), respectively, in the paper (phase 1 = host phase).
 * @param x = r/R	r = The distance from the flat interface, R = The radius of an imaginary sphere.
 * @param pt_in_phase1 = True if the particle is in the particle phase.
 * @return Probability that a test particle first hit \partial\Omega_1 without hitting \partial\Omega_2.*/
	inline double p1_2D(double x, bool pt_in_phase1){
		if(this->alpha > 0.0){
			/* General case */
			if(x < 1e-10){
				return 1.0/(1.0+alpha);
			}
			else{
				if (pt_in_phase1){
					return 1.0/(1.0+alpha) * (1.0 + 4.0*alpha / pi * atan(x));
				}
				else{
					return 1.0/(1.0+alpha) * (1.0 - 4.0 / pi * atan(x));
				}
			}
		}
		else{
			/* if one phase is insulating, 
			 zero probability to cross the interface. */
			if(pt_in_phase1)
				return 1.0;
			else
				return 0.0;
		}
	}

	/** The mean time that a test particle in phase 1 (\theta = \pi/2) or phase 2 (\theta = 3\pi/2) hits \partial \Omega in 2D for the first time.
	 * @param r	The distance from the flat interface, 
	 * @param R The radius of an imaginary sphere* 
	 * @param pt_in_phase1 = True if the particle is in the particle phase.
	 * @return The mean time that a test particle first hit the boundary of the imaginary disk.*/
	inline double tau_s_2D(double r, double R, bool pt_in_phase1){
		double x = r/R;
		if(this->alpha > 0.0){
			double tau_1 = R*R/(4.0*sigma_host);
			if(x < 1e-10){
				return tau_1 * 2.0 / (1.0 + alpha);
			}
			else{
				double rR_Rr = x + 1.0/x;
				if (pt_in_phase1){
					double corr = (1.0 - alpha * x*x +(alpha-1.0)/pi * (atan(x)*rR_Rr*rR_Rr + (x-1.0/x) ) );
					return tau_1* 2.0/(1.0+alpha)* corr;
				}
				else{
					double corr = (1.0 - x*x/alpha -(alpha-1.0)/pi/alpha * (atan(x)*rR_Rr*rR_Rr + (x-1.0/x) ) );
					return tau_1* 2.0/(1.0+alpha)* corr;
				}
			}
		}
		else{
			/* one phase is insulating. */
			if ( x < 1e-10){
				return R*R/ 2.0 /(this->sigma_host + this->sigma_particle);
			}
			else{
				double rR_Rr = x + 1.0/x;
				double corr = (1.0 -1.0/pi * (atan(x)*rR_Rr*rR_Rr + (x-1.0/x) ) );

				return R*R/ 2.0 /(this->sigma_host + this->sigma_particle)* corr;
			}
		}
	}


	/* 2D version for packings of polygonal cells. */


	void PerformSimulation2D_test(const Dispersion & config, size_t num_steps, Trajectory & trajectory, RandomGenerator & rng);

	void PerformSimulation2D_time(const Dispersion & config, double tau_max, Trajectory & trajectory, RandomGenerator & rng);

	void PerformSimulation2D_parallel(std::function<const Dispersion(size_t)> GetConfiguration, size_t num_configs, double time_max, std::vector<GeometryVector> & results, const std::string & log_filename = "");

	void PerformSimulation2D_TimeEvolution(const Dispersion & config, double time_max, double time_interval, std::vector<GeometryVector> & results, const std::string & log_filename = "");

	/* 2D version for disk(sphere) packings. */

	void PerformSimulation2D_test(const SpherePacking & config, size_t num_steps, Trajectory & trajectory, RandomGenerator & rng);

	void PerformSimulation2D_time(const SpherePacking & config, double tau_max, Trajectory & trajectory, RandomGenerator & rng);

	void PerformSimulation2D_parallel(std::function<const SpherePacking(size_t)> GetConfiguration, size_t num_configs, double time_max, std::vector<GeometryVector> & results, const std::string & log_filename = "");

	void PerformSimulation2D_TimeEvolution(const SpherePacking & config, double time_max, double time_interval, std::vector<GeometryVector> & results, const std::string & log_filename = "");

protected:
	void JumpOneStep_2D(const Dispersion & config, Trajectory & trajectory, RandomGenerator & rng, std::vector<double> & debug_info, bool debug = false);

	void JumpOneStep_2D(const SpherePacking & config, Trajectory & trajectory, RandomGenerator & rng, std::vector<double> & debug_info, bool debug = false);
};

#endif