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

/**
 * @brief Find the vertices from two configurations that are separated by lareger than rmin.
 * 
 * @param ifile 
 * @param ofile 
 * @return int 
 */
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


/**
 * @brief Find the vertices from two configurations that are separated by lareger than rmin.
 * Different from ConfigDiff(...), it displaces the compared configuration by the mean displacements
 * 
 * @param ifile 
 * @param ofile 
 * @return int 
 */
int ConfigDiff2(std::istream & ifile, std::ostream & ofile){
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

	double rmean = rmin; //0.1 * pow(comp.PeriodicVolume()/comp.NumParticle(),1./comp.GetDimension());
	std::vector<bool> IdxOfComp_common(comp.NumParticle(), false);
	std::vector<GeometryVector> displacements(comp.NumParticle(), GeometryVector(0,0,0,0));
	comp.PrepareIterateThroughNeighbors(0.1);

	ofile << "Find displacements \n";
	#pragma omp parallel for 
	for (int i=0; i<ref.NumParticle(); i++){
		size_t j;
		bool found = false;
		double dR = -1.0;
		comp.IterateThroughNeighbors(
			ref.GetRelativeCoordinates(i), 0.05,
			[&j, &found](const GeometryVector &x0, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void {
				found = true;
				j = prt_index;
			}, &found		
		);
		if (found){
			GeometryVector displacement =ref.GetRelativeCoordinates(i) - comp.GetRelativeCoordinates(j);
			ref.RelativeCoordToMinimumImage(displacement);
			dR = ref.RelativeCoord2CartesianCoord(displacement).Modulus2();
			IdxOfComp_common[j] = true;
			displacements.at(j) = GeometryVector(dR);
			//# pragma atomic
			//MeanSquaredDisplacement_w_rmin += dR;
		}
		else{
			GeometryVector temp = ref.GetCartesianCoordinates(i);
			temp.Dimension = ::MaxDimension;
			temp.x[::MaxDimension-1] = i;
			//#pragma critical
			//{
			//	OnlyInRef.emplace_back(temp);
			//}
		}
	}
	ofile << "\n";
	
	/* Determine a translation vector (=The most probable one)*/
	GeometryVector x0;
	{	
		struct  comparison
		{
			bool operator() (const GeometryVector& lhs, const GeometryVector& rhs) const
			{
				GeometryVector x = lhs - rhs;
				return   (x.x[0] < 0.0) && (x.x[1] < 0.0);}
			/* data */
		};

		/* Make a histogram */		
		std::map<GeometryVector, int, comparison> distribution;
		for (size_t i =0; i < IdxOfComp_common.size(); i++){
			if ( IdxOfComp_common[i]) {
				distribution[displacements[i]] ++;
			}
		}
		/* Take the most probable displacement to be the translation vector. */
		int highest_count = 0;
		for (auto & x : distribution){
			if (x.second > highest_count){
				highest_count = x.second;
				x0 = GeometryVector(x.first);
			}
		}
	}

	ofile << "Translation vector = " << std::setprecision(14)<< x0 << std::endl;

	// Translate the comparator and then find vertices that are not moved.
	IdxOfComp_common = std::vector<bool>(comp.NumParticle(), false);
	displacements = std::vector<GeometryVector>(comp.NumParticle(), GeometryVector(0,0,0,0));
	for(int i=0; i< comp.NumParticle(); i++){
		comp.MoveParticle(i, comp.GetRelativeCoordinates(i) + x0);
	}
	#pragma omp parallel for reduction(+:MeanSquaredDisplacement_w_rmin)
	for (int i=0; i<ref.NumParticle(); i++){
		size_t j;
		bool found = false;
		double dR = -1.0;
		comp.IterateThroughNeighbors(
			ref.GetRelativeCoordinates(i), rmean,
			[&j, &found](const GeometryVector &x0, const GeometryVector &L_shift, const signed long *PeriodicShift, size_t prt_index) ->void {
				found = true;
				j = prt_index;
			}, &found		
		);
		if (found){
			GeometryVector displacement =ref.GetRelativeCoordinates(i) - comp.GetRelativeCoordinates(j);
			ref.RelativeCoordToMinimumImage(displacement);
			dR = ref.RelativeCoord2CartesianCoord(displacement).Modulus2();
			IdxOfComp_common[j] = true;
			displacements.at(j) = GeometryVector(dR);
		}
		else{
			GeometryVector temp = ref.GetCartesianCoordinates(i);
			temp.Dimension = ::MaxDimension;
			temp.x[::MaxDimension-1] = i;
		}
		MeanSquaredDisplacement_w_rmin += dR;
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
	sprintf(tempstring, "Translation of the comparator from the reference pattern = (%0.13e, %0.13e)", x0.x[0], x0.x[1]);
	fields.emplace_back("coordinates\t\t\t       particle index");
	WriteFunction(OnlyInRef, (savename+"_OnlyInRef").c_str(), fields);
	WriteFunction(OnlyInComp, (savename+"_OnlyInComp").c_str(), fields);

	return 0;
}


#endif // __CONFIGDIFF_H__