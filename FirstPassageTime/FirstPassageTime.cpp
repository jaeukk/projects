/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: April 2021 */

/** \file FirstPassageTime.cpp 
 * \brief Function implementations for the FirstPassageTime simulations. */

#include "FirstPassageTime.h"

const double delta1 = 0.0001;
const double delta2 = 0.01;
const double delta3 = 0.05; // factor for disk packings

const double FACTOR = 1.0 - delta1 * 1e-3;
/* ----------------------------------------------
		Trajectory class
 *----------------------------------------------- */
Trajectory::Trajectory(const Trajectory & src){
	this->dimension = src.dimension;
	this->sigma_host = src.sigma_host; this->sigma_particle = src.sigma_particle;
	this->R2_host = src.R2_host;	this->R2_particle = src.R2_particle;
	this->R2_s = src.R2_s; tau_s = src.tau_s;
	this->tau_total = src.tau_total;

	this->pos_curr = GeometryVector(src.pos_curr);
	this->delta_pos = GeometryVector(src.delta_pos);
}

void Trajectory::ReadLog(std::istream & ifile){
	std::string phrase("#d");
	size_t dim = 0;

	std::string buffer; std::istringstream iss;

	while (dim == 0 && std::getline(ifile, buffer)){
		if(buffer[0] == '#'){
			dim = (size_t) std::stoi(buffer.substr(phrase.size()));	
		}
	}

	if (dim == this->dimension){
		// Read pos_curr
		std::getline(ifile, buffer);
		this->pos_curr.SetDimension(dim);
		iss.str(buffer);
		iss.clear();
		for (int j=0; j<dim; j++)
			iss >> this->pos_curr.x[j];
		
//		if (ifile.eof() == false)
//			vertices.push_back(temp);

		// Read delta_pos
		std::getline(ifile, buffer);
		this->delta_pos.SetDimension(dim);
		iss.str(buffer);
		iss.clear();
		for (int j=0; j<dim; j++)
			iss >> this->delta_pos.x[j];


		std::getline(ifile, buffer);
		iss.str(buffer);
		iss.clear();
		double temp1 = 0.0, temp2 = 0.0;
		// Read R2_host
		iss >> temp1;
		this->Update_R2_host(temp1);

		// Read R2_particle
		iss >> temp1;
		this->Update_R2_particle(temp1);

		// Read R2_s, tau_s
		iss >> temp1; iss >> temp2;
		this->Update_R2_s(temp1, temp2);
		
		// compare time
		iss >> temp1;
		double err_time = temp1 - this->GetTimeTotal();
		this->tau_total += err_time;
	}
	else{
		std::cerr<<"Trajectory::ReadLog(...)\twrong dimension!\n";
	}
}

std::ostream & operator << (std::ostream & os, const Trajectory & t) 
{
//	os << "[Trajectory]\n";
	os << "#d\t" << t.dimension<<"\n";
	os << t.pos_curr;
	os << t.delta_pos;

	os << t.R2_host << "\t" << t.R2_particle <<"\t";
	os << t.R2_s << "\t" << t.tau_s << "\t";
	os << t.tau_total<<"\n";
	return os;
}


/* ------------------------------------------------
	FirstPassageTimeParameters
 * ------------------------------------------------*/

bool FirstPassageTimeParameters::ContinueTrajectory(Trajectory & traj) 
{
	if (! ((traj.sigma_particle == this->sigma_particle)
		&& (traj.sigma_host == this->sigma_host) && (traj.dimension == this->dimension))){
		
		/* Initialize the trajectory */

		traj = Trajectory();
		traj.sigma_particle = this->sigma_particle;
		traj.sigma_host = this->sigma_host;
		traj.dimension = this->dimension;
		traj.pos_curr = GeometryVector(static_cast<int>(traj.dimension));
		traj.delta_pos = GeometryVector(static_cast<int>(traj.dimension));

		return false;
	}
	else{
		return true;
	}
}


void FirstPassageTimeParameters::WriteLogs(const std::string & name, const std::vector<Trajectory> & trajs) 
{
	std::fstream output((name+".trajlog").c_str(), std::fstream::out);

	output <<"##TrajectoryLog##\n";
	output <<"#dimension\n"<<this->dimension<<"\n";
	output <<"#sigma_host\n"<<this->sigma_host<<"\n";
	output <<"#sigma_particle\n"<<this->sigma_particle<<"\n";
	
	size_t num = std::min(this->num_test_prts, trajs.size());
	output <<"#num_test_particles\n"<< num <<"\n";
	output.precision(10);
	output << std::scientific;
	for (size_t i=0; i<num; i++){
		output << "[Trajectory]\n";
		output << trajs[i];
	}

	output.close();
};

bool FirstPassageTimeParameters::ReadLogs(const std::string & name, std::vector<Trajectory> & trajs){

	std::fstream ifile ((name+".trajlog").c_str(), std::fstream::in);

	if (ifile.good() == false)
	{
//		std::cerr << "Error in FirstPassageTimeParameters::ReadLogs : log file doesn't exist. generate a new one.\n" << name << "\n";
		return false;
	}
//	std::cout<<" Continue from the past trajectories\n";

	std::string buffer;
	std::getline(ifile, buffer);
	if (strcmp(buffer.c_str(), "##TrajectoryLog##") == 0){
		DimensionType dimension_read = 0;
		double sigma1 = 0.0, sigma2 = 0.0;
		size_t num_prt_read = 0;

		ifile >> buffer;
		if(strcmp(buffer.c_str(), "#dimension") == 0){
			ifile >> dimension_read;
		}

		ifile >> buffer;
		if(strcmp(buffer.c_str(), "#sigma_host") == 0){
			ifile >> sigma1;
		}

		ifile >> buffer;
		if(strcmp(buffer.c_str(), "#sigma_particle") == 0){
			ifile >> sigma2;
		}

		if (! ((dimension_read == this->dimension) 
			&& (sigma1 == this->sigma_host) 
			&& (sigma2 == this->sigma_particle) ) ){
			std::cerr << "Error in FirstPassageTimeParameters::ReadLogs : log information is different from the current object !\n" << name << "\n";
			return false;
		}

		ifile >> buffer;
		if(strcmp(buffer.c_str(), "#num_test_particles") == 0){
			ifile >> num_prt_read;
		}

		if (num_prt_read < this->num_test_prts){
			std::cerr << "Error in FirstPassageTimeParameters::ReadLogs : the number of trajectories in log is smaller than that in the current object !\n" << name << "\n";
			return false;
		}
		else{
			for (size_t i=0; i< this->num_test_prts; i++){
				ifile >> buffer;
				if(strcmp(buffer.c_str(), "[Trajectory]") == 0){
					trajs[i].ReadLog(ifile);
				}
			}
		}
		return true;
	}
	else{
		return false;
	}
}

void FirstPassageTimeParameters::JumpOneStep_2D(const Dispersion & config, Trajectory & traj, RandomGenerator & rng, std::vector<double> & debug_info, bool doDebug) 
{
	/* Assume that trajectory is correctly initialized. */

	/* Debug Information  
	 *	debug_info[0] = x-coordinate of the position of the test particle. (CartesianCoordinates)
	 *	debug_info[1] = y-coordinate of the position of the test particle. (CartesianCoordinates)
	 *	debug_info[2] = jump radius
	 *	debug_info[3] = current phase 
	 *	debug_info[4] = p1 (probability to go to the host phase, i.e., phase 1)	 */
	debug_info.clear();
	debug_info.resize(this->dimension + 3);

	/* Initial position of a test particle. */
	GeometryVector x(traj.pos_curr);
	
	bool in_host_phase = (config.IsInsideParticle(x) == -1); 

	std::vector<GeometryVector> list;
	GeometryVector x0, dump, displacement;
	config.FindNearestNeighbors(x, list, x0, dump);
	double r = list[0].x[0]; /* the distance to the closest interface */
	double R; 	/* jump radius */
	double p_host, tau_s = 0.0;

	/* -----------------
		Case 1: a test particle is far way from the interface.
	   -----------------*/
	if(r > delta1){
		/* In bulk, */
		//R = r;
		R = FACTOR * r;	// FACTOR <~1 is used to prevent errors due to significant figures.
		displacement = R * RandomUnitVector(this->dimension, rng);
		x = x + config.CartesianCoord2RelativeCoord(displacement);

		if(in_host_phase){

			#ifdef SAFETY_CHECK
			if(!(config.IsInsideParticle(x) == -1)){
				std::cout << ": error!\n";
				GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
				std::cout << "x_prev = " << x_prev;
				std::cout << "x_new = " << x;
			}
			#endif

			traj.Update_R2_host(R*R);
//			traj.R2_host += R*R;
			p_host = 1;
		}
		else{

			#ifdef SAFETY_CHECK
			if((config.IsInsideParticle(x) == -1)){
				std::cout << ": error!!\n";
				GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
				std::cout << "x_prev = " << x_prev;
				std::cout << "x_new = " << x;
			}
			#endif
			
//			traj.R2_particle += R*R;
			traj.Update_R2_particle(R*R);
			p_host = 0;
		}

	}
	else{
	/* -----------------
		Case 2: a test particle is near the interface.
	   -----------------*/		
		GeometryVector normal2host_Cartesian;
		config.FindNearestNeighbors(x0, list, dump, normal2host_Cartesian);
		R = FACTOR * list[1].x[0];	
		/*R is the second nearest distnace from x0 . 
		 * FACTOR <~1 is used to prevent errors due to significant figures.		*/
		bool go2host;


		if (R > delta1){
			/* -------------
				Case 2-1: the interface is near flat.
			   ------------- */

			p_host = this->p1_2D(r/R, in_host_phase); 
			
			double prob = rng.RandomDouble();
			go2host = prob < p_host;
			displacement = R*RandomUnitVector(this->dimension, rng);
			GeometryVector x_trial_Relative = x0 + config.CartesianCoord2RelativeCoord(displacement);

			bool is_trial_in_host = (config.IsInsideParticle(x_trial_Relative)==-1);

			if (go2host != is_trial_in_host){
				// flip the displacement in the normal vector direction!!
				displacement = displacement - 2.0* normal2host_Cartesian.Dot(displacement) * normal2host_Cartesian;
				x_trial_Relative = x0 + config.CartesianCoord2RelativeCoord(displacement);
			}

			tau_s = this->tau_s_2D(r, R,  in_host_phase);

			in_host_phase = go2host;
			x = x_trial_Relative;

			#ifdef SAFETY_CHECK
			if (go2host){
				//Check whether a test particle went to the host phase
				if(!(config.IsInsideParticle(x) == -1)){
					std::cout << " error in near interface case (interface-particle phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}					
			}
			else{
				if (config.IsInsideParticle(x) == -1){
					std::cout << " error in near interface case (interface-host phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}
			}
			#endif
		}
		else{
			/* -------------------
				Case 2-2: Interface of more than two particles are close -> 				The interface has a sharp corner.
			   ------------------- */
			/* In this case, we take advantage of the following approximate formula: 
			
			   p_1 ~ A_1 \sigma_1 / \sum A_i \sigma_i
			
			\tau_s ~ \tau_1 (\sum A_i \sigma_1 / \sum A_i \sigma_i) 			*/
	
			R = delta2;
			r = list[1].x[0];
			size_t num_samples = 60;
			
			/* A simple estimation of the areas in phase 1 and 2. */
			std::vector<GeometryVector> x_particle, x_host;
			x_particle.reserve(num_samples);
			x_host.reserve(num_samples);

			for (int j = 0; j < num_samples; j++){
				double theta = 2.0*pi * (double) j /(double) num_samples + 1e-6;

				GeometryVector x_on_shell = x0 + config.CartesianCoord2RelativeCoord( R*GeometryVector(cos(theta), sin(theta)) );

				if (config.IsInsideParticle(x_on_shell)==-1){
					x_host.push_back(x_on_shell);
				}
				else{
					x_particle.push_back(x_on_shell);
				}
			} 

			double fraction_particle = (double)x_particle.size() / (double) num_samples;
			double fraction_host = (double)x_host.size() / (double) num_samples;

			p_host = fraction_host * this->sigma_host / (fraction_particle * this->sigma_particle + fraction_host * this->sigma_host);

			/* Choose the position of the test particle. */
			double prob = rng.RandomDouble();
			go2host = (prob < p_host);
			if(go2host){
				x = x_host[(size_t)floor(x_host.size() * rng.RandomDouble())];
				in_host_phase = true;
			}
			else{
				x = x_particle[(size_t)floor(x_particle.size() * rng.RandomDouble())];
				in_host_phase = false;
			}
			displacement = config.RelativeCoord2CartesianCoord(x-x0);


			/* Compute the mean hitting time */
			tau_s =  R*R/4.0 / (fraction_particle * this->sigma_particle + fraction_host * this->sigma_host);
			// tau_s *= (1.0 - 0.5*(this->alpha + 1.0) * (r/R)*(r/R) + 0.5*(this->alpha - 1.0) * (r/R)*(r/R)* (r/R)*(r/R));

			#ifdef SAFETY_CHECK
			if (go2host){
				if (!(config.IsInsideParticle(x) == -1)){
					std::cout << " error in near interface case (corner-particle phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}
			}
			else{
				if (config.IsInsideParticle(x) == -1){
					std::cout << " error in near interface case (corner-host phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}
			}
			#endif
		}


//		traj.R2_s += R*R;
//		traj.tau_s += tau_s;
		traj.Update_R2_s(R*R, tau_s);
	}

	if(doDebug){
		GeometryVector x_center(traj.pos_curr);
		if (tau_s > 0.0){
			x_center = x0;
		}
		config.RelativeCoordToMinimumImage(x_center);
		for (size_t i =0; i<this->dimension; i++){
			x_center.x[i] -= std::floor(x_center.x[i]);
			x_center.x[i] -= std::floor(x_center.x[i]);
		}
		x_center = config.RelativeCoord2CartesianCoord(x_center);

		size_t i = 0;
		for (; i<this->dimension; i++){
			debug_info[i] = x_center.x[i];
		}
		debug_info[i++] = R;
		debug_info[i++] = (in_host_phase)? 1:2;
		debug_info[i++] = p_host;
	}

	// Update the total displacement
	traj.delta_pos = traj.delta_pos + displacement;
	traj.pos_curr = x;
}

//For a disk packing
void FirstPassageTimeParameters::JumpOneStep_2D(const SpherePacking & config, Trajectory & traj, RandomGenerator & rng, std::vector<double> & debug_info, bool doDebug) 
{
	/* Assume that trajectory is correctly initialized. */

	/* Debug Information  
	 *	debug_info[0] = x-coordinate of the position of the test particle. (CartesianCoordinates)
	 *	debug_info[1] = y-coordinate of the position of the test particle. (CartesianCoordinates)
	 *	debug_info[2] = jump radius
	 *	debug_info[3] = current phase 
	 *	debug_info[4] = p1 (probability to go to the host phase, i.e., phase 1)	 */
	debug_info.clear();
	debug_info.resize(this->dimension + 3 + 2);

	/* Initial position of a test particle. */
	GeometryVector x(traj.pos_curr);
	
	bool in_host_phase = (config.IsInsideParticle(x) == -1); 

	std::vector<GeometryVector> list;
	GeometryVector x0, dump, displacement;
	config.FindNearestNeighbors(x, list, x0, dump);
	double r = list[0].x[0]; /* the distance to the closest interface */
	double R; 	/* jump radius */
	double p_host, tau_s = 0.0;

	/* -----------------
		Case 1: a test particle is far way from the interface.
	   -----------------*/
	if(r > delta1){
		/* In bulk, */
		//R = r;
		R = FACTOR * r;	// FACTOR <~1 is used to prevent errors due to significant figures.
		displacement = R * RandomUnitVector(this->dimension, rng);
		x = x + config.CartesianCoord2RelativeCoord(displacement);

		if(in_host_phase){

			#ifdef SAFETY_CHECK
			if(!(config.IsInsideParticle(x) == -1)){
				std::cout << ": error!\n";
				GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
				std::cout << "x_prev = " << x_prev;
				std::cout << "x_new = " << x;
			}
			#endif

			traj.Update_R2_host(R*R);
//			traj.R2_host += R*R;
			p_host = 1;
		}
		else{

			#ifdef SAFETY_CHECK
			if((config.IsInsideParticle(x) == -1)){
				std::cout << ": error!!\n";
				GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
				std::cout << "x_prev = " << x_prev;
				std::cout << "x_new = " << x;
			}
			#endif
			
//			traj.R2_particle += R*R;
			traj.Update_R2_particle(R*R);
			p_host = 0;
		}

	}
	else{
	/* -----------------
		Case 2: a test particle is near the interface.
	   -----------------*/		
		GeometryVector normal2host_Cartesian;
		config.FindNearestNeighbors(x0, list, dump, normal2host_Cartesian);
		R = std::min(
			FACTOR * list[1].x[0], //1)
			delta3 * config.GetCharacteristics(list[0].x[2])	//2)
		);	
		/* R is the smallest one among 
		 1) the second nearest distnace from x0, 
		 FACTOR <~1 is used to prevent errors due to significant figures.		
		 2) delta3*a
		 */
		bool go2host;


		if (R > delta1){
			/* -------------
				Case 2-1: the interface is near flat.
			   ------------- */

			p_host = this->p1_2D(r/R, in_host_phase); 
			
			double prob = rng.RandomDouble();
			go2host = prob < p_host;
			displacement = R*RandomUnitVector(this->dimension, rng);
			GeometryVector x_trial_Relative = x0 + config.CartesianCoord2RelativeCoord(displacement);

			bool is_trial_in_host = (config.IsInsideParticle(x_trial_Relative)==-1);

			if (go2host != is_trial_in_host){
				// flip the displacement in the normal vector direction!!
				displacement = displacement - 2.0* normal2host_Cartesian.Dot(displacement) * normal2host_Cartesian;
				x_trial_Relative = x0 + config.CartesianCoord2RelativeCoord(displacement);
			}
			dump = normal2host_Cartesian;
			tau_s = this->tau_s_2D(r, R,  in_host_phase);

			in_host_phase = go2host;
			x = x_trial_Relative;

			#ifdef SAFETY_CHECK
			if (go2host){
				//Check whether a test particle went to the host phase
				if(!(config.IsInsideParticle(x) == -1)){
					std::cout << " error in near interface case (interface-particle phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}					
			}
			else{
				if (config.IsInsideParticle(x) == -1){
					std::cout << " error in near interface case (interface-host phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}
			}
			#endif
		}
		else{
			/* -------------------
				Case 2-2: Interface of more than two particles are close -> 				The interface has a sharp corner.
			   ------------------- */
			/* In this case, we take advantage of the following approximate formula: 
			
			   p_1 ~ A_1 \sigma_1 / \sum A_i \sigma_i
			
			\tau_s ~ \tau_1 (\sum A_i \sigma_1 / \sum A_i \sigma_i) 			*/
	
			R = delta2;
			r = list[1].x[0];
			size_t num_samples = 60;
			
			/* A simple estimation of the areas in phase 1 and 2. */
			std::vector<GeometryVector> x_particle, x_host;
			x_particle.reserve(num_samples);
			x_host.reserve(num_samples);

			for (int j = 0; j < num_samples; j++){
				double theta = 2.0*pi * (double) j /(double) num_samples + 1e-6;

				GeometryVector x_on_shell = x0 + config.CartesianCoord2RelativeCoord( R*GeometryVector(cos(theta), sin(theta)) );

				if (config.IsInsideParticle(x_on_shell)==-1){
					x_host.push_back(x_on_shell);
				}
				else{
					x_particle.push_back(x_on_shell);
				}
			} 

			double fraction_particle = (double)x_particle.size() / (double) num_samples;
			double fraction_host = (double)x_host.size() / (double) num_samples;

			p_host = fraction_host * this->sigma_host / (fraction_particle * this->sigma_particle + fraction_host * this->sigma_host);

			/* Choose the position of the test particle. */
			double prob = rng.RandomDouble();
			go2host = (prob < p_host);
			if(go2host){
				x = x_host[(size_t)floor(x_host.size() * rng.RandomDouble())];
				in_host_phase = true;
			}
			else{
				x = x_particle[(size_t)floor(x_particle.size() * rng.RandomDouble())];
				in_host_phase = false;
			}
			displacement = config.RelativeCoord2CartesianCoord(x-x0);


			/* Compute the mean hitting time */
			tau_s =  R*R/4.0 / (fraction_particle * this->sigma_particle + fraction_host * this->sigma_host);
			// tau_s *= (1.0 - 0.5*(this->alpha + 1.0) * (r/R)*(r/R) + 0.5*(this->alpha - 1.0) * (r/R)*(r/R)* (r/R)*(r/R));

			#ifdef SAFETY_CHECK
			if (go2host){
				if (!(config.IsInsideParticle(x) == -1)){
					std::cout << " error in near interface case (corner-particle phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}
			}
			else{
				if (config.IsInsideParticle(x) == -1){
					std::cout << " error in near interface case (corner-host phase)\n";
					GeometryVector x_prev = x - config.CartesianCoord2RelativeCoord(displacement);
					std::cout << "x_prev = " << x_prev;
					std::cout << "x_new = " << x;
				}
			}
			#endif
		}


//		traj.R2_s += R*R;
//		traj.tau_s += tau_s;
		traj.Update_R2_s(R*R, tau_s);
	}

	if(doDebug){
		GeometryVector x_center(traj.pos_curr);
		if (tau_s > 0.0){
			x_center = x0;
		}
		config.RelativeCoordToMinimumImage(x_center);
		for (size_t i =0; i<this->dimension; i++){
			x_center.x[i] -= std::floor(x_center.x[i]);
			x_center.x[i] -= std::floor(x_center.x[i]);
		}
		x_center = config.RelativeCoord2CartesianCoord(x_center);

		size_t i = 0;
		for (; i<this->dimension; i++){
			debug_info[i] = x_center.x[i];
		}
		debug_info[i++] = R;
		debug_info[i++] = (in_host_phase)? 1:2;
		debug_info[i++] = p_host;
		debug_info[i++] = dump.x[0];
		debug_info[i++] = dump.x[1];
	}

	// Update the total displacement
	traj.delta_pos = traj.delta_pos + displacement;
	traj.pos_curr = x;
}


void FirstPassageTimeParameters::PerformSimulation2D_test(const Dispersion & config, size_t num_steps, Trajectory & traj, RandomGenerator & rng){
	assert(config.GetDimension() == 2);

	if (! this->ContinueTrajectory(traj)){
		for (size_t i=0; i<traj.dimension; i++)
			traj.pos_curr.x[i] = rng.RandomDouble();		
	}
	/* Initial position of a test particle. */
	GeometryVector x(traj.pos_curr);
	
	/* Check whether the test particle is stuck in the insulating phase. */
	bool stucked = false;
	if (this->sigma_particle * this->sigma_host == 0.0){
		bool in_host_phase = (config.IsInsideParticle(x) == -1); 

		if(in_host_phase){
			stucked = this->sigma_host == 0.0;
		}
		else{
			stucked = this->sigma_particle == 0.0;
		}
	}

	if(stucked){
		std::cout << "a test particle is stucked in an insulating phase!!\n";
	}

	if(! stucked){
		std::vector<double> debug_info;
		std::fstream ofile("/home/jaeukk/Projects/Networks/transport_properties/dielectric/FirstPassageTime/log1.txt", std::fstream::out);
		
		for(size_t i=0; i<num_steps; i++){
			this->JumpOneStep_2D(config, traj, rng, debug_info, true);

			for (int j=0; j<debug_info.size(); j++)
				ofile << debug_info[j]<<"\t";
			ofile << "\n";
		}

		ofile.close();

	}
}


void FirstPassageTimeParameters::PerformSimulation2D_time(const Dispersion & config, double tau_max, Trajectory & traj, RandomGenerator & rng) 
{
	assert(config.GetDimension() == 2);

	if (! this->ContinueTrajectory(traj)){
		for (size_t i=0; i<traj.dimension; i++)
			traj.pos_curr.x[i] = rng.RandomDouble();		
	}
	/* Initial position of a test particle. */
	GeometryVector x(traj.pos_curr);
	
	/* Check whether the test particle is stuck in the insulating phase. */
	bool stucked = false;
	if (this->sigma_particle * this->sigma_host == 0.0){
		bool in_host_phase = (config.IsInsideParticle(x) == -1); 

		if(in_host_phase){
			stucked = this->sigma_host == 0.0;
			if(stucked){
				traj.Update_R2_host(0.0, tau_max);
			}
		}
		else{
			stucked = this->sigma_particle == 0.0;
			if(stucked){
				traj.Update_R2_particle(0.0, tau_max);
			}
		}
	}

	if(!stucked){
		std::vector<double> debug_info;
		
		double tau_target = traj.GetTimeTotal() + tau_max;
		
		while(traj.GetTimeTotal() < tau_target){
			this->JumpOneStep_2D(config, traj, rng, debug_info, false);
		}
	}

	
}


void FirstPassageTimeParameters::PerformSimulation2D_parallel(std::function<const Dispersion(size_t)> GetConfiguration, size_t num_configs, double time_max, std::vector<GeometryVector> & results, const std::string & log_filename) 
{
	omp_set_num_threads(this->num_threads);
	results.clear();
	results.resize(num_configs, GeometryVector(0.0,0.0,0.0,0.0));

	if (! (strcmp(log_filename.c_str(), "")==0)){
		std::cout << "Load log files\n";
	}
	else{
		std::cout << "No log files\n";
	}
	

	std::cout<< "Computing sigma_e via first-passage-time simulation\n";
	progress_display pd(num_configs);
	for (size_t i=0; i<num_configs; i++){
		char logname[300]={};
		sprintf(logname, "%s-%d", log_filename.c_str(), i);

		std::vector<Dispersion> c(this->num_threads, GetConfiguration(i));
		std::vector<RandomGenerator> RNG;
		{
			double Rmax = c[0].GetMaxRadius();
			for (int j=0; j<this->num_threads; j++){
				c[j].SetCellSize(Rmax);
				RNG.emplace_back(j);
			}
		}
		std::vector<Trajectory> trajs (this->num_test_prts, Trajectory(this->dimension, this->sigma_host, this->sigma_particle));
		/* Load trajectories from a log file, if it is prescribed. */
		bool success_load = false;
		if (! (strcmp(log_filename.c_str(), "")==0)){
			//std::fstream ifile ((log_filename+".trajlog").c_str(), std::fstream::in);
			success_load = this->ReadLogs(logname, trajs);
		}

		if (!success_load){
#pragma omp parallel for schedule(guided)
			for (int j=0; j<this->num_test_prts; j++){
				size_t tid = omp_get_thread_num();
				for (size_t k=0; k<this->dimension; k++){
					trajs[j].pos_curr.x[k] = RNG[tid].RandomDouble();
				}
			}
		}

		std::vector<GeometryVector> data_config (this->num_test_prts, GeometryVector(0.0,0.0,0.0,0.0));
#pragma omp parallel for schedule(guided)
		for (int j = 0; j<this->num_test_prts; j++){
			size_t tid = omp_get_thread_num();
			this->PerformSimulation2D_time(c[tid], time_max, trajs[j], RNG[tid]);

			data_config[j].x[0] = trajs[j].GetTimeTotal();
			data_config[j].x[1] = trajs[j].GetX2();
		}


		/* treat data_config */
		double mean_tau = 0.0, mean_X2 = 0.0;
		double ste_tau = 0.0, ste_X2 = 0.0;

#pragma omp parallel for reduction(+ : mean_tau, mean_X2)
		for(size_t j=0; j<data_config.size(); j++){
			mean_tau += data_config[j].x[0];
			mean_X2 += data_config[j].x[1];
		}
		mean_tau /= (double) data_config.size();
		mean_X2 /= (double) data_config.size();

#pragma omp parallel for reduction(+ : ste_tau, ste_X2)
		for(size_t j=0; j<data_config.size(); j++){
			double temp = (data_config[j].x[0] - mean_tau);
			ste_tau += temp*temp;
			temp = (data_config[j].x[1] - mean_X2);
			ste_X2 += temp*temp;
		}
		ste_tau /= (double)this->num_test_prts;
		ste_tau = sqrt(ste_tau/((double)this->num_test_prts - 1.0));
		ste_X2 /= (double)this->num_test_prts;
		ste_X2 = sqrt(ste_X2/((double)this->num_test_prts - 1.0));

		results[i].x[0] = mean_tau;
		results[i].x[1] = mean_X2;
		results[i].x[2] = ste_tau;
		results[i].x[3] = ste_X2;

		/* Write a log file, if the name is prescribed. */
		if (! (strcmp(log_filename.c_str(), "")==0)){
			this->WriteLogs(logname, trajs);
		}

		pd++;
	}
}

void FirstPassageTimeParameters::PerformSimulation2D_TimeEvolution(const Dispersion & config, double time_max, double time_interval, std::vector<GeometryVector> & results, const std::string & log_filename) 
{
	omp_set_num_threads(this->num_threads);

	size_t num_time_slices = (size_t)floor(time_max / time_interval);
	results.clear();
	results.resize(num_time_slices, GeometryVector(0.0,0.0,0.0,0.0));


	std::vector<Dispersion> c(this->num_threads, config);
	std::vector<RandomGenerator> RNG;
	{
		double Rmax = c[0].GetMaxRadius();
		for (int j=0; j<this->num_threads; j++){
			c[j].SetCellSize(Rmax);
			RNG.emplace_back(j);
		}
	}
	std::vector<Trajectory> trajs (this->num_test_prts, Trajectory(this->dimension, this->sigma_host, this->sigma_particle));
	/* Load trajectories from a log file, if it is prescribed. */
	bool success_load = false;
	if (! (strcmp(log_filename.c_str(), "")==0)){
		//std::fstream ifile ((log_filename+".trajlog").c_str(), std::fstream::in);
		success_load = this->ReadLogs(log_filename, trajs);
	}

	if (!success_load){
#pragma omp parallel for schedule(guided)
		for (int j=0; j<this->num_test_prts; j++){
			size_t tid = omp_get_thread_num();
			for (size_t k=0; k<this->dimension; k++){
				trajs[j].pos_curr.x[k] = RNG[tid].RandomDouble();
			}
		}
	}

	std::cout<< "Computing sigma_e via first-passage-time simulation\n";
	progress_display pd(num_time_slices);

	for (size_t i=0; i<num_time_slices; i++){

		std::vector<GeometryVector> data_config (this->num_test_prts, GeometryVector(0.0,0.0,0.0,0.0));
#pragma omp parallel for schedule(guided)
		for (int j = 0; j<this->num_test_prts; j++){
			size_t tid = omp_get_thread_num();
			this->PerformSimulation2D_time(c[tid], time_interval, trajs[j], RNG[tid]);

			data_config[j].x[0] = trajs[j].GetTimeTotal();
			data_config[j].x[1] = trajs[j].GetX2();
		}

		/* treat data_config */
		double mean_tau = 0.0, mean_X2 = 0.0;
		double ste_tau = 0.0, ste_X2 = 0.0;

#pragma omp parallel for reduction(+ : mean_tau, mean_X2)
		for(size_t j=0; j<data_config.size(); j++){
			mean_tau += data_config[j].x[0];
			mean_X2 += data_config[j].x[1];
		}
		mean_tau /= (double) data_config.size();
		mean_X2 /= (double) data_config.size();

#pragma omp parallel for reduction(+ : ste_tau, ste_X2)
		for(size_t j=0; j<data_config.size(); j++){
			double temp = (data_config[j].x[0] - mean_tau);
			ste_tau += temp*temp;
			temp = (data_config[j].x[1] - mean_X2);
			ste_X2 += temp*temp;
		}
		ste_tau /= (double)this->num_test_prts;
		ste_tau = sqrt(ste_tau/((double)this->num_test_prts - 1.0));
		ste_X2 /= (double)this->num_test_prts;
		ste_X2 = sqrt(ste_X2/((double)this->num_test_prts - 1.0));

		results[i].x[0] = mean_tau;
		results[i].x[1] = mean_X2;
		results[i].x[2] = ste_tau;
		results[i].x[3] = ste_X2;

		pd++;

		/* Write a log file, if the name is prescribed. */
		if (! (strcmp(log_filename.c_str(), "")==0)){
			//std::fstream ofile ((log_filename+".trajlog").c_str(), std::fstream::out);
			this->WriteLogs(log_filename, trajs);
		}
	}
}

/* For disk packings... */
void FirstPassageTimeParameters::PerformSimulation2D_test(const SpherePacking & config, size_t num_steps, Trajectory & traj, RandomGenerator & rng){
	assert(config.GetDimension() == 2);

	if (! this->ContinueTrajectory(traj)){
		for (size_t i=0; i<traj.dimension; i++)
			traj.pos_curr.x[i] = rng.RandomDouble();		
	}
	/* Initial position of a test particle. */
	GeometryVector x(traj.pos_curr);
	
	/* Check whether the test particle is stuck in the insulating phase. */
	bool stucked = false;
	if (this->sigma_particle * this->sigma_host == 0.0){
		bool in_host_phase = (config.IsInsideParticle(x) == -1); 

		if(in_host_phase){
			stucked = this->sigma_host == 0.0;
		}
		else{
			stucked = this->sigma_particle == 0.0;
		}
	}

	if(stucked){
		std::cout << "a test particle is stucked in an insulating phase!!\n";
	}

	if(! stucked){
		std::vector<double> debug_info;
		std::fstream ofile("/home/jaeukk/Projects/Networks/transport_properties/dielectric/FirstPassageTime/log-disks.txt", std::fstream::out);
		
		for(size_t i=0; i<num_steps; i++){
			this->JumpOneStep_2D(config, traj, rng, debug_info, true);

			for (int j=0; j<debug_info.size(); j++)
				ofile << debug_info[j]<<"\t";
			ofile << "\n";
		}

		ofile.close();

	}
}

void FirstPassageTimeParameters::PerformSimulation2D_time(const SpherePacking & config, double tau_max, Trajectory & traj, RandomGenerator & rng){
	assert(config.GetDimension() == 2);

	if (! this->ContinueTrajectory(traj)){
		for (size_t i=0; i<traj.dimension; i++)
			traj.pos_curr.x[i] = rng.RandomDouble();		
	}
	/* Initial position of a test particle. */
	GeometryVector x(traj.pos_curr);
	
	/* Check whether the test particle is stuck in the insulating phase. */
	bool stucked = false;
	if (this->sigma_particle * this->sigma_host == 0.0){
		bool in_host_phase = (config.IsInsideParticle(x) == -1); 

		if(in_host_phase){
			stucked = this->sigma_host == 0.0;
			if(stucked){
				traj.Update_R2_host(0.0, tau_max);
			}
		}
		else{
			stucked = this->sigma_particle == 0.0;
			if(stucked){
				traj.Update_R2_particle(0.0, tau_max);
			}
		}
	}

	if(!stucked){
		std::vector<double> debug_info;
		
		double tau_target = traj.GetTimeTotal() + tau_max;
		
		while(traj.GetTimeTotal() < tau_target){
			this->JumpOneStep_2D(config, traj, rng, debug_info, false);
		}
	}
}

void FirstPassageTimeParameters::PerformSimulation2D_parallel(std::function<const SpherePacking(size_t)> GetConfiguration, size_t num_configs, double time_max, std::vector<GeometryVector> & results, const std::string & log_filename){
	omp_set_num_threads(this->num_threads);
	results.clear();
	results.resize(num_configs, GeometryVector(0.0,0.0,0.0,0.0));

	if (! (strcmp(log_filename.c_str(), "")==0)){
		std::cout << "Load log files\n";
	}
	else{
		std::cout << "No log files\n";
	}
	

	std::cout<< "Computing sigma_e via first-passage-time simulation\n";
	progress_display pd(num_configs);
	for (size_t i=0; i<num_configs; i++){
		char logname[300]={};
		sprintf(logname, "%s-%d", log_filename.c_str(), i);

		std::vector<SpherePacking> c(this->num_threads, GetConfiguration(i));
		std::vector<RandomGenerator> RNG;
		{
			double Rmax = c[0].GetMaxRadius();
			for (int j=0; j<this->num_threads; j++){
				c[j].SetCellSize(Rmax);
				RNG.emplace_back(j);
			}
		}
		std::vector<Trajectory> trajs (this->num_test_prts, Trajectory(this->dimension, this->sigma_host, this->sigma_particle));
		/* Load trajectories from a log file, if it is prescribed. */
		bool success_load = false;
		if (! (strcmp(log_filename.c_str(), "")==0)){
			//std::fstream ifile ((log_filename+".trajlog").c_str(), std::fstream::in);
			success_load = this->ReadLogs(logname, trajs);
		}

		if (!success_load){
#pragma omp parallel for schedule(guided)
			for (int j=0; j<this->num_test_prts; j++){
				size_t tid = omp_get_thread_num();
				for (size_t k=0; k<this->dimension; k++){
					trajs[j].pos_curr.x[k] = RNG[tid].RandomDouble();
				}
			}
		}

		std::vector<GeometryVector> data_config (this->num_test_prts, GeometryVector(0.0,0.0,0.0,0.0));
#pragma omp parallel for schedule(guided)
		for (int j = 0; j<this->num_test_prts; j++){
			size_t tid = omp_get_thread_num();
			this->PerformSimulation2D_time(c[tid], time_max, trajs[j], RNG[tid]);

			data_config[j].x[0] = trajs[j].GetTimeTotal();
			data_config[j].x[1] = trajs[j].GetX2();
		}


		/* treat data_config */
		double mean_tau = 0.0, mean_X2 = 0.0;
		double ste_tau = 0.0, ste_X2 = 0.0;

#pragma omp parallel for reduction(+ : mean_tau, mean_X2)
		for(size_t j=0; j<data_config.size(); j++){
			mean_tau += data_config[j].x[0];
			mean_X2 += data_config[j].x[1];
		}
		mean_tau /= (double) data_config.size();
		mean_X2 /= (double) data_config.size();

#pragma omp parallel for reduction(+ : ste_tau, ste_X2)
		for(size_t j=0; j<data_config.size(); j++){
			double temp = (data_config[j].x[0] - mean_tau);
			ste_tau += temp*temp;
			temp = (data_config[j].x[1] - mean_X2);
			ste_X2 += temp*temp;
		}
		ste_tau /= (double)this->num_test_prts;
		ste_tau = sqrt(ste_tau/((double)this->num_test_prts - 1.0));
		ste_X2 /= (double)this->num_test_prts;
		ste_X2 = sqrt(ste_X2/((double)this->num_test_prts - 1.0));

		results[i].x[0] = mean_tau;
		results[i].x[1] = mean_X2;
		results[i].x[2] = ste_tau;
		results[i].x[3] = ste_X2;

		/* Write a log file, if the name is prescribed. */
		if (! (strcmp(log_filename.c_str(), "")==0)){
			this->WriteLogs(logname, trajs);
		}

		pd++;
	}
}

void FirstPassageTimeParameters::PerformSimulation2D_TimeEvolution(const SpherePacking & config, double time_max, double time_interval, std::vector<GeometryVector> & results, const std::string & log_filename){
	omp_set_num_threads(this->num_threads);

	size_t num_time_slices = (size_t)floor(time_max / time_interval);
	results.clear();
	results.resize(num_time_slices, GeometryVector(0.0,0.0,0.0,0.0));


	std::vector<SpherePacking> c(this->num_threads, config);
	std::vector<RandomGenerator> RNG;
	{
		double Rmax = c[0].GetMaxRadius();
		for (int j=0; j<this->num_threads; j++){
			c[j].SetCellSize(Rmax);
			RNG.emplace_back(j);
		}
	}
	std::vector<Trajectory> trajs (this->num_test_prts, Trajectory(this->dimension, this->sigma_host, this->sigma_particle));
	/* Load trajectories from a log file, if it is prescribed. */
	bool success_load = false;
	if (! (strcmp(log_filename.c_str(), "")==0)){
		//std::fstream ifile ((log_filename+".trajlog").c_str(), std::fstream::in);
		success_load = this->ReadLogs(log_filename, trajs);
	}

	if (!success_load){
#pragma omp parallel for schedule(guided)
		for (int j=0; j<this->num_test_prts; j++){
			size_t tid = omp_get_thread_num();
			for (size_t k=0; k<this->dimension; k++){
				trajs[j].pos_curr.x[k] = RNG[tid].RandomDouble();
			}
		}
	}

	std::cout<< "Computing sigma_e via first-passage-time simulation\n";
	progress_display pd(num_time_slices);

	for (size_t i=0; i<num_time_slices; i++){

		std::vector<GeometryVector> data_config (this->num_test_prts, GeometryVector(0.0,0.0,0.0,0.0));
#pragma omp parallel for schedule(guided)
		for (int j = 0; j<this->num_test_prts; j++){
			size_t tid = omp_get_thread_num();
			this->PerformSimulation2D_time(c[tid], time_interval, trajs[j], RNG[tid]);

			data_config[j].x[0] = trajs[j].GetTimeTotal();
			data_config[j].x[1] = trajs[j].GetX2();
		}

		/* treat data_config */
		double mean_tau = 0.0, mean_X2 = 0.0;
		double ste_tau = 0.0, ste_X2 = 0.0;

#pragma omp parallel for reduction(+ : mean_tau, mean_X2)
		for(size_t j=0; j<data_config.size(); j++){
			mean_tau += data_config[j].x[0];
			mean_X2 += data_config[j].x[1];
		}
		mean_tau /= (double) data_config.size();
		mean_X2 /= (double) data_config.size();

#pragma omp parallel for reduction(+ : ste_tau, ste_X2)
		for(size_t j=0; j<data_config.size(); j++){
			double temp = (data_config[j].x[0] - mean_tau);
			ste_tau += temp*temp;
			temp = (data_config[j].x[1] - mean_X2);
			ste_X2 += temp*temp;
		}
		ste_tau /= (double)this->num_test_prts;
		ste_tau = sqrt(ste_tau/((double)this->num_test_prts - 1.0));
		ste_X2 /= (double)this->num_test_prts;
		ste_X2 = sqrt(ste_X2/((double)this->num_test_prts - 1.0));

		results[i].x[0] = mean_tau;
		results[i].x[1] = mean_X2;
		results[i].x[2] = ste_tau;
		results[i].x[3] = ste_X2;

		pd++;

		/* Write a log file, if the name is prescribed. */
		if (! (strcmp(log_filename.c_str(), "")==0)){
			//std::fstream ofile ((log_filename+".trajlog").c_str(), std::fstream::out);
			this->WriteLogs(log_filename, trajs);
		}
	}
}

