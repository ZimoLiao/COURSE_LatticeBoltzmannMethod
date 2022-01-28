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
	const int host = 0;
	int rank, size;

	// geometry (grid)
	int nx, ny, tg[5] = { 0 };
	int ni, nj;

	// simulation
	int nforce = 5;
	int nstep = 0, nwrite = 0;

	// model
	bool is_srt = false, is_mrt = false;
	double tau = 0.8;

	// boundary counter
	int nlb_total = 0;


	/* internal functions */
	void Ignore(ifstream& fin, int l);

public:
	/* data */
	LatticeMoment lm;
	LatticePopulation lp;


	/* constructor */
	LatticeSolver();

};


#endif // !LATTICESOLVER_H_