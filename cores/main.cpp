#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

//-------Core part--------
#include "etc.h"
#include "GeometryVector.h"
#include "PeriodicCellList.h"
#include "RandomGenerator.h"

//------Analysis part-----
//#include "StructureFactor.h"
//#include "PairCorrelation.h"

//void WriteToContourf(const std::vector<GeometryVector> & vec_k, const std::vector<GeometryVector> &s_k, const std::string &prefix, int length);


#include <set>

int main() {

	char tempstring[1000];
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	char name[200] = {};
	char dir[] = "./results/n1/";
	

	//declared in "RandomGenerator.h"
	//random number generator
	RandomGenerator rng(0);
	for(int i=0; i< 10; i++)
		ofile << rng.RandomDouble() <<"\t";	//generate a real random number from 0 to 1
	
	//declared in "GeometryVector.h"
	//implement a vector quantity in Euclidean space up to d=4
	{
		GeometryVector v1((int) 1);	//v1 = a zero vector in 1D
		GeometryVector v2(1.0);		//v2 = (1.0) in 1D
		GeometryVector v3(rng.RandomDouble(), rng.RandomDouble());	//v3 = a random vector in 2D
		GeometryVector v4(0.5, 0.1, 0.3);	//v4 = (0.5, 0.1, 0.3) in 3D
		
		//print
		ofile << "v1="<<v1 << "v2=" << v2 << "v3="<< v3 << "v4="<< v4;
		
		//operations
		//+
		ofile << "v1+v2=" << v1+v2;
		//multiplication
		ofile << "0.1* v3=" << 0.1*v3;
		//Inner product
		ofile << "Inner prodcut of v4 and (0, 2, 3) = "<< v4.Dot(GeometryVector(0.0, 2.0, 3.0))<<std::endl;				
	}
	
	// PeriodicCellList: declared in "PeriodicCellList.h"
	// implement the periodic simulation box in d-dimensional Euclidean space d<=4
	// neighbor cell list is implemented.
	{
		//Configuration class is a special type of PeriodicCellList for a point configuration
		int dimension = 3;
		std::vector<GeometryVector> basis_vector;
		{	basis_vector.emplace_back(2.0, 0.0, 0.0);
			basis_vector.emplace_back(0.0, 1.0, 0.0);
			basis_vector.emplace_back(0.0, 0.0, 2.0);}
		double CellSize = 0.2; 
		Configuration con(dimension, &basis_vector[0], CellSize);	
		// con = a periodic simulation box in d=3, 
		//	 basis_vector describes the basis vectors
		//	 CellSize is an approximate size of an individual cell
		
		// Insert 10 particles in "con"
		// Particle coordinates are interpreted as relative coordinates in terms of basis vectors.
		// Whenever a particle is inserted, its is recorded in three forms:
		// (i) Cartesian Coordinates
		// (ii) Relative Coordinates
		// (iii) Cell in which the particle is located.
		int particle_number = 100;
		for(int i=0; i<particle_number; i++){
			GeometryVector new_relative_coordinate(rng.RandomDouble(), rng.RandomDouble(), rng.RandomDouble());
			con.Insert("a", new_relative_coordinate);	//"a" is the name of particle
		}
		
		//interfaces
		ofile 	<< "\n\n\n\n"
			<< "-----------------------------\n"
			<< "        Configuration        \n"
			<< "-----------------------------\n"
			<< "Dimension = " << con.GetDimension()<<"\n"
			<< "basis vector[0] = "<< con.GetBasisVector(0)
			<< "basis vector[1] = "<< con.GetBasisVector(1)
			<< "basis vector[2] = "<< con.GetBasisVector(2)<<"\n\n"
			<< "There are "<< con.NumParticle() << " in this object\n"
			<< "Print cartesian coordinates of 10 particles among them\n\n";
		for(int i=0; i< 10; i++)
			ofile << i <<": "<<con.GetCartesianCoordinates(i);

		ofile <<"con is saved in \"con.txt\"\n";
		WriteConfiguration(con, "con");
	}	
		
	/*  Compute the isotropic structure factor */
	//Lambda function to load configurations
	int number_configuration = 1;
	double Kmax_circular = 10.0, Kmax_linear = 100.0, dK = 2.0;  
	std::string filename("con");
	auto GetConfiguration = [&filename](size_t i)->Configuration{	
		Configuration temp;	// an empty Configuration object
		ReadConfiguration(temp, filename);	
		//there are multiple configurations with different indices one can use the following expression
		//ReadConfiguration(temp, (filename+std::to_string(i)).c_str())	;	
		return temp;
	};	
	
	std::vector<GeometryVector> result_Sk;
	/*----------------------------------------
	   IsotropicStructureFactor(...) is declared in "StructureFactor.h"
	   Compute isotropic structure factor in the following steps:
	   Find all wavevectors whose magnitude is smaller than 'Kmax_circular'
	// Then, obtain all wavevectors that satisfy  |(integer)*k| < 'Kmax_linear' from wavevectors k in the previous step
	// Compute structure factor at all wavevectors that we find.
	// Average structure factor according to magnitude of wavevectors, where bin width is 'dK'
	// dK must be larger than 2*pi/(box length)

	   IsotropicSpectralDensity(...) is the counterpart of the spectral density for "sphere packings"	
	*/
	// IsotropicStructureFactor(GetConfiguration, number_configuration, Kmax_circular, Kmax_linear,result_Sk, dK );
	// ofile << "S(k) is saved in Sk_con.txt\n";
	// WriteFunction(result_Sk, "Sk_con");

				
	/*  Compute isotropic pair correlation function */
	// double R_max = 1.0;
	// std::vector<GeometryVector> result_g2;
	// IsotropicTwoPairCorrelation(GetConfiguration, number_configuration, R_max, result_g2);	
	// ofile << "g_2(r) is saved in g2_con.txt\n";
	// WriteFunction(result_g2, "g2_con");

	return 0;
}


// void WriteToContourf(const std::vector<GeometryVector> & vec_k, const std::vector<GeometryVector> &s_k, const std::string &prefix, int length) {

// 	if (vec_k.size() == s_k.size()) {
// 		if (vec_k[0].Dimension == 2) {
// 			std::string file = prefix + "X_.txt";
// 			std::ofstream X_(file);
// 			file = prefix + "Y_.txt";
// 			std::ofstream Y_(file);
// 			file = prefix + "S_k.txt";
// 			std::ofstream S(file);

// 			for (int i = 0; i < vec_k.size(); i++) {
// 				X_ << vec_k[i].x[0] << "\t";
// 				Y_ << vec_k[i].x[1] << "\t";
// 				S << s_k[i].x[0] << "\t";
// 				if ((i + 1) % length == 0) {
// 					X_ << "\n";
// 					Y_ << "\n";
// 					S << "\n";
// 				}
// 			}
// 			X_.close();	Y_.close(); S.close();
// 		}
// 		else if (vec_k[0].Dimension == 3) {
// 			std::string file = prefix + "X_.txt";
// 			std::ofstream X_(file);
// 			file = prefix + "Y_.txt";
// 			std::ofstream Y_(file);
// 			file = prefix + "Z_.txt";
// 			std::ofstream Z_(file);
// 			file = prefix + "S_k.txt";
// 			std::ofstream S(file);

// 			for (int i = 0; i < vec_k.size(); i++) {
// 				X_ << vec_k[i].x[0] << "\t";
// 				Y_ << vec_k[i].x[1] << "\t";
// 				Z_ << vec_k[i].x[2] << "\t";
// 				S << s_k[i].x[0] << "\t";
// 				if ((i + 1) % length == 0) {
// 					X_ << "\n";
// 					Y_ << "\n";
// 					Z_ << "\n";
// 					S << "\n";
// 				}
// 			}
// 			X_.close();	Y_.close(); Z_.close();  S.close();
// 		}
// 	}
// 	else {
// 		std::cout << " wrong input \n";
// 	}



// }
