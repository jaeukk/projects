/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: October 2024 */

/** \file main.cpp 
 * \brief A simple CLI for performing some calculations for QC configurations*/


#include "ConfigDiff.h"

/* header fields in $(cores) */

size_t Verbosity = 2;

int main(int argc, char ** argv){
	char tempstring[1000] = {};
	std::istream &ifile = std::cin;
	std::ostream &ofile = std::cout;

	if (argc > 1 && strcmp (argv[1], "ConfigDiff") == 0){
		ConfigDiff(ifile, ofile);
	}	

	return 0;	
}
