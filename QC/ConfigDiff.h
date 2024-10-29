#ifndef __CONFIGDIFF_H__
#define __CONFIGDIFF_H__

/* header fields in $(cores) */
#include <PeriodicCellList.h>
#include <GeometryVector.h>
#include <etc.h>


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

bool comparison (size_t i,size_t j) { return (i<=j); }

int ConfigDiff(std::istream & ifile, std::ostream & ofile){
	char tempstring[1000] = {};
	Configuration ref, comp;
	double rmin = 0;
	double MeanSquaredDisplacement_w_rmin = 0.0;
	size_t num_threads;
	std::vector<GeometryVector> OnlyInRef, OnlyInComp;
	std::string savename;

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

		ofile << "_______________________\n";
		ofile << "Compared configuration\n";
		ofile << "-----------------------\n";
		ofile << "ConfigPack Name:\n";
		ifile >> tempstring;
		ofile << tempstring << std::endl;
		ConfigurationPack Comp(tempstring);
		ofile << "Index:\n";
		ifile >> idx;
		comp = Comp.GetConfig(idx);
	}

	/* Set the minimum radius to distinguish particle positions*/
	{
		ofile << "R_min = ";
		ifile >> rmin;
		ofile << rmin << std::endl;
		ofile << "Number of threads = ";
		ifile >> num_threads;
		ofile << num_threads << std::endl;
		omp_set_num_threads(num_threads);
		ofile << "save name = ";
		ifile >> savename;
		ofile << savename << std::endl;
	}

	std::vector<bool> IdxOfComp_common(comp.NumParticle(), false);
	comp.PrepareIterateThroughNeighbors(rmin);
	
	ofile << "start comparison\n";
	#pragma omp parallel for 
	for (int i=0; i<ref.NumParticle(); i++){
		size_t j;
		bool found = false;
		double dR = -1.0;
		comp.IterateThroughNeighbors(
			ref.GetRelativeCoordinates(i), rmin,
			[&j, &found](const GeometryVector &x0, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void {
				j = prt_index;
				found = true;
			}		
		);
		if (found){
			GeometryVector displacement =ref.GetRelativeCoordinates(i) - comp.GetRelativeCoordinates(j);
			ref.RelativeCoordToMinimumImage(displacement);
			dR = ref.RelativeCoord2CartesianCoord(displacement).Modulus2();
			IdxOfComp_common[j] = true;
			
			# pragma atomic
			MeanSquaredDisplacement_w_rmin += dR;
		}
		else{
			GeometryVector temp = ref.GetCartesianCoordinates(i);
			temp.Dimension = ::MaxDimension;
			temp.x[::MaxDimension-1] = i;
			#pragma critical
			{
				OnlyInRef.emplace_back(temp);
			}
		}
	}

	ofile << "Comparisons are done\n";

	size_t j = 0;
	double count = 0;
	for(size_t i=0; i<comp.NumParticle(); i++){
		if(! IdxOfComp_common[i]){
			GeometryVector temp = comp.GetCartesianCoordinates(i);
			temp.Dimension = ::MaxDimension;
			temp.x[::MaxDimension-1] = i;
			OnlyInComp.emplace_back(temp);
		}
		else{
			count ++ ;
		}
	}

	/* Save !! */
	if (IdxOfComp_common.size() > 0){
		MeanSquaredDisplacement_w_rmin /= count;
	}
	else{
		MeanSquaredDisplacement_w_rmin = 0.0L;
	}

	sprintf(tempstring, "Mean Squared Displacement = %0.10e", MeanSquaredDisplacement_w_rmin);
	ofile << tempstring << std::endl;
	std::vector<std::string> fields;
	fields.emplace_back(std::string(tempstring)+"\n");
	fields.emplace_back("coordinates\t\t\t       particle index");
	WriteFunction(OnlyInRef, (savename+"_OnlyInRef").c_str(), fields);
	WriteFunction(OnlyInComp, (savename+"_OnlyInComp").c_str(), fields);

	return 0;
}




#endif // __CONFIGDIFF_H__