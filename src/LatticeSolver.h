#ifndef LATTICESOLVER_H_
#define LATTICESOLVER_H_

#include<mpi.h>
#include<iostream>
#include<fstream>
#include<string>

#include"LatticeBound.h"
#include"LatticeMoment.h"
#include"LatticePopulation.h"

using namespace std;

class LatticeSolver
{
	/* parameters */
	int rank, size;

	// geometry (grid)
	int nx, ny, tg[5] = { 0,1,1,1,1 };
	int ni, nj;

	// simulation
	int nforce = 5;
	int nstep = 0;

public:
	/* data */
	LatticeMoment lm;
	LatticePopulation lp;


	/* constructor */
	LatticeSolver();

};


#endif // !LATTICESOLVER_H_