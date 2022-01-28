#ifndef LATTICESOLVER_H_
#define LATTICESOLVER_H_

#include<mpi.h>
#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include<time.h>

#include"LatticeBound.h"
#include"LatticeMoment.h"
#include"LatticePopulation.h"
#include"ParticleSolver.h"

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
	int nforce = 1;
	int nstep = 0, nwrite = 0;

	// model
	bool is_srt = false, is_mrt = false;
	double tau = 0.8;

	// boundary counter
	int nlb_total = 0;


	/* particle solver */
	ParticleSolver psolver;


	/* internal functions */
	void Ignore(ifstream& fin, int l);

public:
	/* data */
	LatticeMoment lm;
	LatticePopulation lp;


	/* monitor (global) */
	int step = 0;
	double diff;

	/* constructor */
	LatticeSolver();

	/* functions */
	void Calculate();

};


#endif // !LATTICESOLVER_H_