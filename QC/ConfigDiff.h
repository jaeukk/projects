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
	bool PrintDebugData = false;
	
	std::string reference, comparator;

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
		reference = std::string(tempstring);

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
		comparator = std::string(tempstring);
	}

	/* Set the minimum radius to distinguish particle positions*/
	{
		// ofile << "R_min = ";
		// ifile >> rmin;
		// ofile << rmin << std::endl;
		ofile << "Number of threads = ";
		ifile >> num_threads;
		ofile << num_threads << std::endl;
		omp_set_num_threads(num_threads);
		ofile << "save name = ";
		ifile >> savename;
		ofile << savename << std::endl;

	}

	//double rmean = 0.05; //rmin; //0.1 * pow(comp.PeriodicVolume()/comp.NumParticle(),1./comp.GetDimension());
	std::vector<bool> IdxOfComp_common(comp.NumParticle(), false);
	std::vector<GeometryVector> displacements(comp.NumParticle(), GeometryVector(0,0,0,0));
	double TypicalLength = std::pow(comp.PeriodicVolume() / comp.NumParticle(), 1.0 / comp.GetDimension())*1.6; 
	comp.PrepareIterateThroughNeighbors(TypicalLength*3.0);

	ofile << "Find displacements \n";
	{
		//#pragma omp parallel for 
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
				IdxOfComp_common[j] = true;
				displacements.at(j) = GeometryVector(displacement); //dR);
			}
		}
		ofile << "\n";
	}
	
	/* Determine a translation vector (=The most probable one)*/
	GeometryVector x0 = GeometryVector(0.0L, 0.0L);
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
		if (PrintDebugData){
			WriteFunction(displacements, (savename+"_Debug").c_str());
		}
	}


	ofile << "Translation vector = " << std::setprecision(14)<< x0 << std::endl;
	
	// Translate the comparator by x0.
	for(int i=0; i< comp.NumParticle(); i++){
		comp.MoveParticle(i, comp.GetRelativeCoordinates(i) + x0);
	}
	// Now, x0 is in cartesian coordinates.
	x0 = comp.RelativeCoord2CartesianCoord(x0);	

	// Find the vertices that are not moved.
	ofile << "Compare the updated comparator to the reference.\n";
	displacements = std::vector<GeometryVector>(comp.NumParticle(), GeometryVector(0,0,0,0));	// from a vertex in reference.
	std::vector<size_t> indices_in_comp(ref.NumParticle());

	#pragma omp parallel for reduction(+:MeanSquaredDisplacement_w_rmin)
	for (int i = 0; i < ref.NumParticle(); i++ ){
		/* Find the closest particle of the comparator for particle i in the reference. */
		std::vector<GeometryVector> neighbors;
		std::vector<size_t> idx_neighbors;
		GeometryVector ref_pt = ref.GetRelativeCoordinates(i);

		double l = TypicalLength;
		while (neighbors.size() < 1)
		{
			neighbors.clear();
			idx_neighbors.clear();

			comp.IterateThroughNeighbors(ref_pt, l, [&neighbors, &idx_neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
				idx_neighbors.push_back(Sourceparticle);
			});
			l *= 1.5;
		}

		double nearest_displacement = 10000.0;
		GeometryVector dx;
		size_t idx_comp;
		for(size_t j = 0; j < neighbors.size(); j++ ){
			double val = neighbors[j].Modulus2();
			if (nearest_displacement > val){
				nearest_displacement = val;
				dx = GeometryVector(neighbors[j]);
				idx_comp = idx_neighbors[j];
			}
		}
		//dx = comp.RelativeCoord2CartesianCoord(dx);
		indices_in_comp.at(i) = idx_comp;	
		displacements.at(i).x[0] = dx.x[0];
		displacements.at(i).x[1] = dx.x[1];
		displacements.at(i).x[3] = (double) i;

		/* Check whether the displacement is smaller than rmin */
		MeanSquaredDisplacement_w_rmin += dx.Modulus2();
	}
	MeanSquaredDisplacement_w_rmin /= ref.NumParticle();

	// sorting displacements in the descending order of magnitudes.
	std::sort(displacements.begin(), displacements.end(), 
		[](const GeometryVector & left, const GeometryVector & right)->bool{return (left.x[0]*left.x[0]+left.x[1]*left.x[1]) > (right.x[0]*right.x[0]+right.x[1]*right.x[1]);} 
	);

	/* Save !! */
	std::vector<size_t> idx_for_test;
	idx_for_test.reserve(ref.NumParticle());
	{
		std::fstream save((savename+".txt").c_str(), std::fstream::out);
		/* Headers */
		sprintf(tempstring, "# Reference: %s.ConfigPack, N=%zu \n", reference.c_str(), ref.NumParticle());
		save << tempstring; 
		sprintf(tempstring, "# Comparator: %s.ConfigPack, N=%zu \n", comparator.c_str(), comp.NumParticle());
		save << tempstring; 
		sprintf(tempstring, "# Mean Squared Displacement = %0.10e \n", MeanSquaredDisplacement_w_rmin);
		save << tempstring;
		sprintf(tempstring, "# Translation = (%0.13e, %0.13e) \n", x0.x[0], x0.x[1]);
		save << tempstring; 
		sprintf(tempstring, "# Number density = %0.6e \n", comp.NumParticle()/comp.PeriodicVolume());
		save << tempstring;
		save << "# id1 = particle index in reference \n";
		save << "# id2 = particle index in comparator \n";
		save << "# x = x-coordinate in displacement \n";
		save << "# y = y-coordinate in displacement \n";
		save << "# id1 \t\t id2 \t\t x \t\t y \n";

		/* Body */
		size_t i;
		for (i = 0; i < displacements.size(); i++){
			double val = displacements[i].x[0]*displacements[i].x[0] + displacements[i].x[1]* displacements[i].x[1];
			size_t id = (size_t) displacements[i].x[3];

			if (val > 1e-20){
				// Print pairs of vertices that were moved within a machine precision.
				sprintf(tempstring, "%zu\t%zu\t%0.10e\t%0.10e\n", id, indices_in_comp[id], displacements[i].x[0], displacements[i].x[1]);
				save << tempstring;

				idx_for_test.emplace_back(indices_in_comp[id]);
			}
		}
		save.close();
	}

	// Sanity check: all elements of idx_for_test should be unique!
	{
		size_t diff = idx_for_test.size();
		std::sort(idx_for_test.begin(), idx_for_test.end()); 
	    auto last = std::unique(idx_for_test.begin(), idx_for_test.end());
	    idx_for_test.erase(last, idx_for_test.end());

		diff -= idx_for_test.size();
		if ( diff > 0 ){
			ofile << "There are some duplicates in indices of comparator:\t" << diff << "\n";
		}
	}

	return 0;
}


#endif // __CONFIGDIFF_H__