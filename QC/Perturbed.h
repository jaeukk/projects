#ifndef __PERTURBED_H__
#define __PERTURBED_H__

/* header fields in $(cores) */
#include <PeriodicCellList.h>
#include <GeometryVector.h>
#include <etc.h>
#include <RandomGenerator.h>


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

//bool comparison (size_t i,size_t j) { return (i<=j); }

/**
 * @brief Generate perturbed lattices.
 * 
 * @param ifile 
 * @param ofile 
 * @return int 
 */
int Perturbed(std::istream & ifile, std::ostream & ofile){
	/* Possible perturbations */
	std::vector<GeometryVector> displacements;
	displacements.emplace_back(-0.6180339888,0.0000000000);	//0
	displacements.emplace_back(-0.5000000000,0.3632712640); //1
	displacements.emplace_back(-0.3090169944,-0.2245139883); //2
	displacements.emplace_back(-0.1909830056,-0.5877852523); //3
	displacements.emplace_back(-0.1180339888,0.3632712640); //4
	displacements.emplace_back(0.1909830056,-0.5877852523); //5
	displacements.emplace_back(0.1909830056,0.5877852523); //6
	displacements.emplace_back(0.3819660112,-0.0000000000); //7
	displacements.emplace_back(0.5000000000,-0.3632712640); //8
	displacements.emplace_back(0.5000000000,0.3632712640); //9


	char tempstring[1000] = {};
	Configuration ref;
	size_t num_perturbed, Nc, seed;
	std::string savename;
	std::vector<GeometryVector> dx;
	RandomGenerator rng;
	/* Two configurations to be compared */
	{
		ofile << "_______________________\n";
		ofile << "Reference configuration\n";
		ofile << "-----------------------\n";
		ofile << "ConfigPack Name:\n";
		ifile >> tempstring;
		ofile << tempstring << std::endl;
		ConfigurationPack Ref(tempstring);
		ofile << "Index:\n";
		size_t idx;
		ifile >> idx;
		ref = Ref.GetConfig(idx);
		for (size_t i = 0; i < displacements.size(); i++){
			dx.push_back(ref.CartesianCoord2RelativeCoord(displacements[i]));
		}
	}

	/* Set the minimum radius to distinguish particle positions*/
	{
		ofile << "Number of displaced vertices = ";
		ifile >> num_perturbed;
		ofile << num_perturbed << std::endl;

		ofile << "Number of configurations = ";
		ifile >> Nc;
		ofile << Nc << std::endl;

		ofile << "random seed = ";
		ifile >> seed;
		ofile << seed << std::endl;
		
		ofile << "save name = ";
		ifile >> savename;
		ofile << savename << std::endl;
	}

	ConfigurationPack save(savename);
	rng.seed(seed);
	size_t N = ref.NumParticle();
	for (size_t i = 0; i < Nc; i++){
		Configuration c (ref);
		for (size_t j=0; j <  num_perturbed; j++){
			size_t id = (size_t)std::floor(N * rng.RandomDouble());
			size_t k = (size_t)std::floor(10 * rng.RandomDouble());
			if (id == N)
				id --;
			if (k == 10)
				k -- ;
			c.MoveParticle(id, c.GetRelativeCoordinates(id) + dx[k]);
		}
		save.AddConfig(c);
	}

	ofile << "end!" << std::endl;

	return 0;
}



#endif // __PERTURBED_H__