#ifndef LATTICEPOPULATION_H_
#define LATTICEPOPULATION_H_

#include<mpi.h>
#include<vector>
#include<iostream>

#include"LatticeBound.h"

using std::vector;
using std::cin;
using std::cout;
using std::endl;

class LatticePopulation
{

	/* data */
	// populations array
	double* data;

	// fluid-boundaries
	int nlb;
	vector<LatticeBound> lb;

public:
	/* constructor & destructor */
	~LatticePopulation();

	/* initialization */
	void InitBoundary(LatticeBound newlb);

};


#endif // !LATTICEPOPULATION_H_