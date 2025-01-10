#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>



#include "../etc.h"
#include "../GeometryVector.h"
#include "../RandomGenerator.h"

#include "../NonsphericalParticles/ParticleShapes.h"

#include "../PeriodicCellList.h"
#include "../NonsphericalParticles/Dispersions.h"

#include "FirstPassageTime.h"
#include <chrono>

void WriteFunction(const std::vector<GeometryVector> & result, std::ostream & ofile, const std::vector<std::string> & fields);
void WriteFunction(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::vector<std::string> & fields);
void ReadFunction(std::vector<GeometryVector> & result, std::istream & ifile, size_t NumColumns);
void ReadFunction(std::vector<GeometryVector> & result, const std::string & InFilePrefix, size_t NumColumns);

void Test_IsInside();
void Test_distance();
void Test_Dispersion();
void Test_p1();
void Test_simulation(double x, double y, size_t);
void HoneycombNetwork(double, double, double, double, double);
void Checkerboard(double, double, double, double);

void Test_interface(double x, double y, size_t);

int main(int argc, char ** argv) {
	std::string tempStr;
	char tempstring[1000];
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	RandomGenerator rng(0);
	int num_threads = 6;
	omp_set_num_threads(num_threads);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

//	Trajectory test;
//	ofile << test;
//	double x = std::atof(argv[1]);
//	double y = std::atof(argv[2]);
//	Test_simulation2();
	if(strcmp(argv[1],"0")==0){
		Checkerboard(std::stof(argv[2]),std::stof(argv[3]),std::stof(argv[4]),std::stof(argv[5]));
	}
	else if(strcmp(argv[1],"1")==0){
		HoneycombNetwork(std::stof(argv[2]),std::stof(argv[3]),std::stof(argv[4]),std::stof(argv[5]), std::stof(argv[6]));
	}
	else if(strcmp(argv[1],"2")==0){
		for (int i=0; i<argc; i++)
			ofile << argv[i] <<std::endl;

		std::string loadname (argv[2]);
		size_t start_idx = std::atoi(argv[3]);
		size_t num_configs = std::atoi(argv[4]);
		std::string savename (argv[5]);
		std::string logname("");
		if (argc == 7){
			logname = std::string(argv[6]);
		}
		
		std::vector<std::string> header;
		header.emplace_back("First-passage-time calculation for 2D packings \n");

		sprintf(tempstring, "Load %s\n", loadname.c_str());
		header.emplace_back(tempstring);
		ofile << tempstring;

		sprintf(tempstring, "\tfrom %d to %d\n", start_idx, start_idx + num_configs);
		header.emplace_back(tempstring);
		ofile << tempstring;

		sprintf(tempstring, "save to %s\n", savename.c_str());
		header.emplace_back(tempstring);
		ofile << tempstring;

		if (argc == 7){
			sprintf(tempstring, "save log files as %s\n", logname.c_str());
			header.emplace_back(tempstring);
			ofile << tempstring;
		}

		double sigma1, sigma2;
		double time_max;
		size_t num_test_prts_per_config, num_threads;


		ofile << "sigma in host phase = ";
		ifile >> sigma1;
		ofile << sigma1 <<"\n";
		ofile << "sigma in particle phase = ";
		ifile >> sigma2;
		ofile << sigma2 << "\n";
		ofile << "target simulation time = ";
		ifile >> time_max;
		ofile << time_max <<"\n";
		ofile << "num. test particles / config = ";
		ifile >> num_test_prts_per_config;
		ofile << num_test_prts_per_config <<"\n";
		ofile << "num. threads = ";
		ifile >> num_threads;
		ofile << num_threads <<"\n";


		auto getc = [&loadname, &start_idx](size_t i)->Dispersion{
			Dispersion c;
			char name[400] = {};
			sprintf(name, loadname.c_str(), start_idx + i);
			ReadDispersion(c, name);
			return c;
		};

		FirstPassageTimeParameters p(sigma1, sigma2, getc(0).GetDimension());
		p.num_test_prts = num_test_prts_per_config;
		p.num_threads = num_threads;

		std::vector<GeometryVector> result;
		p.PerformSimulation2D_parallel(getc, num_configs, time_max, result, logname);

		
		header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
		header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
		header.emplace_back("num. random walkers / configs= " + std::to_string(p.num_test_prts) + "\n");
		header.emplace_back("<time>\t<X2>\terr. <time>\terr. <X2>");
		WriteFunction(result, savename, header);
	}
	else if(strcmp(argv[1],"3")==0){
		int num = 500;
		double x = atof(argv[2]);
		double y = atof(argv[3]);
		SpherePacking c;
		ReadConfiguration(c, "./config");
		ofile <<"read\n";
		double count_total_trials = (double) num, count2host = 0.0;
		double sigma1 = 1.0, sigma2 = 0.0;
		ofile << std::setprecision(8);
		ofile << std::scientific;
		// for (double X = 0.3-0.0001; X < 0.6 ; X += 0.05){
		// 	GeometryVector v (X, 0.0), x0_Relative, n2host;
		// 	std::vector<GeometryVector> list;
		// 	c.FindNearestNeighbors(v, list, x0_Relative, n2host);
		
		// 	ofile << "v = \n";
		// 	ofile << v;
		// 	ofile << "list = \n";
		// 	for (auto u = list.begin(); u!=list.end(); u++){
		// 		ofile << * u;
		// 	}
		// 	ofile << "x0 = \n";
		// 	ofile << x0_Relative << "\n";
		// 	ofile << "n2 host = \n";
		// 	ofile << n2host;

		// 	ofile << "\n\n\n";
		// }

		
		FirstPassageTimeParameters param(sigma1, sigma2, 2);
		Trajectory one;
		RandomGenerator rng(0);
		one.sigma_host = sigma1;	one.sigma_particle = sigma2;
		one.dimension = 2;
		one.delta_pos = GeometryVector(0.0,0.0);
		one.pos_curr = GeometryVector(x,y);
		std::cout << "initial position = " << one.pos_curr;
		std::vector<GeometryVector> debug;
		for(size_t i=0; i<1; i++){
			param.PerformSimulation2D_test(c, num, one, rng);
			//if(c.IsInsideParticle(one.pos_curr) == -1){
			//	count2host ++;
			//}

			//one.pos_curr = GeometryVector(x,y);
			//std::cout<<one.pos_curr;
			//param.JumpOneStep_2D(c, one, rng, debug, true);

		}
		ofile << one;
		
		std::cout << one.GetTimeTotal() << "\t" << one.GetX2()<<"\n";
		std::cout<<"p_1 = "<< count2host / count_total_trials <<"\n";
	}
	else if(strcmp(argv[1],"4")==0){
		for (int i=0; i<argc; i++)
			ofile << argv[i] <<std::endl;

		SpherePacking c;
		if(strcmp(argv[2],"square")==0){
			ReadConfiguration(c, "./square");
		}
		else if(strcmp(argv[2],"triangular")==0){
			ReadConfiguration(c, "./triangular");
		}

		std::string name(argv[2]);
		double phi2 = atof(argv[3]);
		double sigma1 = atof(argv[4]);	// sigma_host
		double sigma2 = atof(argv[5]);	// sigma_particle
		double time_max = atof(argv[6]);
		double time_interval = atof(argv[7]);
		std::string logname("");
		if (argc == 9){
			logname=std::string(argv[8]);
		}

		double phi_curr = c.PackingFraction();
		double factor = sqrt(phi2/phi_curr);
		c.Rescale_prtOnly(factor);
		std::cout << "particles are rescaled: packing fraction = "<< c.PackingFraction()<<"\n";
		
		FirstPassageTimeParameters param(sigma1, sigma2, 2);
		param.num_test_prts = 10; //00;
		param.num_threads = 4;

		std::vector<GeometryVector> result;
		param.PerformSimulation2D_TimeEvolution(c, time_max, time_interval, result, logname);

		char filename[400]={};
		std::vector<std::string> header;

		if(strcmp(argv[2],"square")==0){
			sprintf(filename, "/tigress/jaeukk/FPT/square/square-%0.2f-%0.1f-%0.1f", phi2, sigma1, sigma2);
			header.emplace_back("First-passage-time calculation for 2D periodic square lattice packing of identical circles\n");
		}
		else if(strcmp(argv[2],"triangular")==0){
			sprintf(filename, "/tigress/jaeukk/FPT/triangular/triangular-%0.2f-%0.1f-%0.1f", phi2, sigma1, sigma2);
			header.emplace_back("First-passage-time calculation for 2D periodic triangular lattice packing of identical circles\n");
		}
				
		header.emplace_back("area fraction of phase 2 = " + std::to_string(phi2) + "\n");
		header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
		header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
		header.emplace_back("num. random walkers = " + std::to_string(param.num_test_prts) + "\n");
		header.emplace_back("mean time\tmean X2\terr. time\terr. X2");

		WriteFunction(result, filename, header);
	}

//	Test_simulation(std::atof(argv[1]),std::atof(argv[2]), std::atoi(argv[3]));

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	ofile << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

//	RandomGenerator RNG;
//	Test_simulation(RNG.RandomDouble(),RNG.RandomDouble(), num_steps);

//	Test_interface(std::atof(argv[1]),std::atof(argv[2]), num_steps);
//	Test_simulation3(std::atof(argv[1]),std::atof(argv[2]), num_steps);

	//Test_p1();
	// std::vector<GeometryVector> test;
	// for (int i=0; i<1000; i++){
	// 	test.push_back(RandomUnitVector(2, rng));
	// }
	// WriteFunction(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test");

	return 0;
}

void Checkerboard(double time_max, double time_interval, double sigma1, double sigma2){
	
	Dispersion c;
	ReadDispersion(c, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/checkerboard/regular2");
	
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	param.num_test_prts = 100; //00;
	param.num_threads = 1;

	char name[400]={};
	sprintf(name, "./checkerboard-%0.2f-%0.2f", sigma1, sigma2);
	std::vector<GeometryVector> result;
	param.PerformSimulation2D_TimeEvolution(c, time_max, time_interval, result, name);

	sprintf(name, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/checkerboard-%0.1f-%0.1f", sigma1, sigma2);
	std::vector<std::string> header;
	header.emplace_back("First-passage-time calculation for 2D periodic checkerboard\n");
	header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
	header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
	header.emplace_back("num. random walkers = " + std::to_string(param.num_test_prts) + "\n");
	header.emplace_back("mean time\tmean X2\terr. time\terr. X2");

	WriteFunction(result, name, header);
}

void HoneycombNetwork(double time_max, double time_interval, double phi2, double sigma1, double sigma2){
	
	Dispersion c;
	ReadDispersion(c, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
	double phi_curr = c.PackingFraction();
	double factor = sqrt(phi2/phi_curr);
	c.Rescale_prtOnly(factor);
	std::cout << "particles are rescaled: packing fraction = "<< c.PackingFraction()<<"\n";
	
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	param.num_test_prts = 1000000; //00;
	param.num_threads = 20;

	std::vector<GeometryVector> result;
	param.PerformSimulation2D_TimeEvolution(c, time_max, time_interval, result);

	char name[400]={};
	sprintf(name, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/honeycomb-%0.2f-%0.1f-%0.1f", phi2, sigma1, sigma2);
	std::vector<std::string> header;
	header.emplace_back("First-passage-time calculation for 2D periodic honeycomb network\n");
	header.emplace_back("area fraction of phase 2 = " + std::to_string(phi2) + "\n");
	header.emplace_back("simga1 = " + std::to_string(sigma1) + "\n");
	header.emplace_back("simga2 = " + std::to_string(sigma2) + "\n");
	header.emplace_back("num. random walkers = " + std::to_string(param.num_test_prts) + "\n");
	header.emplace_back("mean time\tmean X2\terr. time\terr. X2");

	WriteFunction(result, name, header);
}



void Test_interface(double x, double y, size_t num){

	double count_total_trials = (double) num, count2host = 0.0;
	double sigma1 = 1.0, sigma2 = 4.0;
	
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	Dispersion test;
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/Z2-tri");
	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/checkerboard/regular2");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/2D_HSF_0.20_100-f_0.5-0");
	Trajectory one;
	RandomGenerator rng(0);
	one.sigma_host = sigma1;	one.sigma_particle = sigma2;
	one.dimension = 2;
	one.pos_curr = GeometryVector(x,y);
	std::cout << "initial position = " << one.pos_curr;
	for(size_t i=0; i<num; i++){
		param.PerformSimulation2D_test(test, 1, one, rng);
		if(test.IsInsideParticle(one.pos_curr) == -1){
			count2host ++;
		}
		one.pos_curr = GeometryVector(x,y);
	}
	std::cout<<"p_1 = "<< count2host / count_total_trials <<"\n";

}



void Test_simulation(double x, double y, size_t num_steps){
	double sigma1 = 1.0, sigma2 = 0.0;
	FirstPassageTimeParameters param(sigma1, sigma2, 2);
	Dispersion test;
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/Z2-tri");
	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/checkerboard/regular2");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/2D_HSF_0.20_100-f_0.5-0");
	Trajectory one;
	RandomGenerator rng(0);
	one.sigma_host = sigma1;	one.sigma_particle = sigma2;
	one.dimension = 2;
	one.pos_curr = GeometryVector(x,y);
	one.delta_pos = GeometryVector(0,0);
	std::cout << "initial position = " << one.pos_curr;

	param.PerformSimulation2D_test(test, num_steps, one, rng);

	char temp[300];
	std::cout<<"R2 total, R2_host, R2_prt, R2_s, tau_s\n";
//	sprintf(temp, "%1.7e, %1.7e, %1.7e, %1.7e, %1.7e, %1.7e\n", one.GetR2Total(), one.R2_host, one.R2_particle, one.R2_s, one.tau_s, one.GetEffectiveConductivity());
	std::cout << temp;
	std::cout << one.R2_host / one.R2_particle <<"\n";
//	if (i==0){
//		sprintf(temp, "data = np.array([[%1.7e, %1.7e, %1.7e, %1.7e, %1.7e, %1.7e,]", one.GetR2Total(), one.R2_host, one.R2_particle, one.R2_s, one.tau_s, one.GetEffectiveConductivity());
//	}
//	else{
//		sprintf(temp, ",\n[%1.7e, %1.7e, %1.7e, %1.7e, %1.7e, %1.7e,]", one.GetR2Total(), one.R2_host, one.R2_particle, one.R2_s, one.tau_s, one.GetEffectiveConductivity());
//	}

//	ofile<<temp;
//	ofile<<"])\n";

	 
}

void Test_p1(){
	FirstPassageTimeParameters param(1.0, 10.0, 2);
	double R = 1.0;
	
	std::vector<GeometryVector> p1_1, p1_2, taus_1, taus_2;

	for (double r = 0.01; r<R; r+=0.01){
		p1_1.emplace_back(r, param.p1_2D(r, true));
		p1_2.emplace_back(r, param.p1_2D(r, false));

		taus_1.emplace_back(r, param.tau_s_2D(r, R, true));
		taus_2.emplace_back(r, param.tau_s_2D(r, R, false));
	}
	WriteFunction(p1_1, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/p1_1");
	WriteFunction(p1_2, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/p1_2");
	WriteFunction(taus_1, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/taus_1");
	WriteFunction(taus_2, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/taus_2");
}

void Test_Dispersion(){
	/* read an example dispersion. */
	Dispersion test;
//	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/honeycomb-0.50");
	ReadDispersion(test, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/poly_dispersions/temp/2D_HSF_0.20_100-f_0.5-0");

	std::cout << * test.GetShape(86);
	std::cout << "\n";
	std::cout << test.GetCartesianCoordinates(86);
	RandomGenerator rng(1);
	std::vector<GeometryVector> inside, outside;
	inside.reserve(1000);
	outside.reserve(1000);

	for (int i=0; i<1000; i++){
		GeometryVector p_test(rng.RandomDouble()/10., rng.RandomDouble()/10.);
		GeometryVector p_Cartesian = test.RelativeCoord2CartesianCoord(p_test);
		p_Cartesian.SetDimension(3);
		p_Cartesian.x[2] = i;
		if(test.IsInsideParticle(p_test)==-1){
			outside.push_back(p_Cartesian);
		}
		else{
			inside.push_back(p_Cartesian);
		}
	}

	WriteFunction(inside, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/pts_inside");
	WriteFunction(outside, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/pts_outside");

	std::vector<GeometryVector> centers, origins, radii;


	for (int i=0; i<30; i++){
		GeometryVector p_test(rng.RandomDouble()/10., rng.RandomDouble()/10.);
		GeometryVector p_Cartesian = test.RelativeCoord2CartesianCoord(p_test);
		p_Cartesian.SetDimension(3);
		p_Cartesian.x[2] = i;

		std::vector<GeometryVector> list;
		GeometryVector x0, n;

		test.FindNearestNeighbors(p_test, list, x0, n);
		
		centers.push_back(p_Cartesian);
		origins.push_back(test.RelativeCoord2CartesianCoord(x0));
		radii.emplace_back(list[0].x[0], list[1].x[0], list[2].x[0]);

	}
	WriteFunction(centers, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test_particles");
	WriteFunction(origins, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test_x0");
	WriteFunction(radii, "/home/jaeukk/Projects/OrderMetrics/transport_properties/dielectric/FirstPassageTime/test_radii");
}

void Test_IsInside(){
	std::vector<GeometryVector> vs;
	char temp_str[300]="";
	vs.emplace_back(-0.5, -0.5);
	vs.emplace_back(0.5, -0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(0.4, 0.3);
	vs.emplace_back(-0.5, 0.5);

	Polygon temp(vs);
	std::fstream ofile("test_isinside.txt", std::fstream::out);



	ofile<< "vs = np.array([";
	for (size_t i=0; i < vs.size(); i++){
		sprintf(temp_str, "[%0.6f,%0.6f],\n", vs[i].x[0],vs[i].x[1]);

		ofile<< temp_str;
	}
	sprintf(temp_str, "[%0.6f,%0.6f]])\n\n", vs[0].x[0],vs[0].x[1]);
	ofile << temp_str;


	int num = 1000, count = 0;
	RandomGenerator rng(312);
	std::vector<GeometryVector> inside, outside;
	inside.reserve(num);
	outside.reserve(num);

	GeometryVector cm = temp.GetCenterOfMass();
	for (int i=0; i <num ; i++){
		GeometryVector ptr(2.0*(rng.RandomDouble()-0.5), 2.0*(rng.RandomDouble()-0.5));

		if(temp.IsPtInside(ptr-cm)){
			inside.emplace_back(ptr.x[0],ptr.x[1],(double)i);
		}
		else{
			outside.emplace_back(ptr.x[0],ptr.x[1],(double)i);
		}
	}

	ofile<<"ptr1=np.array([";
	for (int i=0; i<inside.size(); i++){
		if (i == 0)
			sprintf(temp_str, "[%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2]);
		else{
			sprintf(temp_str, ",\n[%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";

	ofile<<"ptr2=np.array([";
	for (int i=0; i<outside.size(); i++){
		if (i == 0)
			sprintf(temp_str, "[%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2]);
		else{
			sprintf(temp_str, ",\n[%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";
	ofile.close();
}


void Test_distance(){
	std::vector<GeometryVector> vs;
	char temp_str[300]="";
	vs.emplace_back(-0.5, -0.5);
	vs.emplace_back(0.5, -0.5);
	vs.emplace_back(0.5, 0.5);
	vs.emplace_back(0.4, 0.7);
	vs.emplace_back(-0.5, 0.5);

	Polygon temp(vs);
	std::fstream ofile("test_distance.txt", std::fstream::out);



	ofile<< "vs = np.array([";
	for (size_t i=0; i < vs.size(); i++){
		sprintf(temp_str, "[%0.6f,%0.6f],\n", vs[i].x[0],vs[i].x[1]);

		ofile<< temp_str;
	}
	sprintf(temp_str, "[%0.6f,%0.6f]])\n\n", vs[0].x[0],vs[0].x[1]);
	ofile << temp_str;


	int num = 1000, count = 0;
	RandomGenerator rng(312);
	std::vector<GeometryVector> inside, outside;
	inside.reserve(num);
	outside.reserve(num);

	GeometryVector cm = temp.GetCenterOfMass();
	for (int i=0; i <num ; i++){
		GeometryVector ptr(2.0*(rng.RandomDouble()-0.5), 2.0*(rng.RandomDouble()-0.5));
		GeometryVector ptr_center = ptr-cm;
		std::vector<GeometryVector> list;
		GeometryVector x0, n;
		double distance = temp.GetDistance2Boundary(ptr_center, list, x0, n);
		if(temp.IsPtInside(ptr_center)){
			inside.emplace_back(ptr.x[0],ptr.x[1], -1.0*distance, (double)i);
		}
		else{
			outside.emplace_back(ptr.x[0],ptr.x[1], distance, (double)i);
		}

		if (i < 10){
			std::cout<< "----\t"<<i <<"\t----\n";
			std::cout<< "x = " << ptr;
			std::cout<< "x0 = "<< x0 + cm;
			std::cout<< "\n";
			for (auto l = list.begin(); l != list.end(); l++)
				std::cout << *l;
			
			
			std::cout<<"\n\n";
		}
	}

	ofile<<"ptr1=np.array([";
	for (int i=0; i<inside.size(); i++){
		if (i == 0)
			sprintf(temp_str, "[%0.6f,%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2], inside[i].x[3]);
		else{
			sprintf(temp_str, ",\n[%0.6f,%0.6f,%0.6f,%0.1f]", inside[i].x[0], inside[i].x[1], inside[i].x[2], inside[i].x[3]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";

	ofile<<"ptr2=np.array([";
	for (int i=0; i<outside.size(); i++){
		if (i == 0)
			sprintf(temp_str, "[%0.6f,%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2],outside[i].x[3]);
		else{
			sprintf(temp_str, ",\n[%0.6f,%0.6f,%0.6f,%0.1f]", outside[i].x[0], outside[i].x[1],outside[i].x[2],outside[i].x[3]);}
		ofile<< temp_str;
	}
	ofile<<"])\n\n";
	ofile.close();
}