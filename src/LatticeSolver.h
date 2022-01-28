#ifndef LATTICESOLVER_H_
#define LATTICESOLVER_H_

#include<mpi.h>
#include<iostream>
#include<fstream>
<<<<<<< HEAD
#include<algorithm>
#include<vector>
#include<string>
#include<iomanip>

#include"LatticeBound.h"
#include"LatticeEntity.h"
#include"LatticeMoment.h"
#include"LatticePopulation.h"

using namespace std;

/* Lattice Boltzmann solver for D2Q9 model
	Properties:

	Author:
		Zimo Liao,
		Department of Modern Mechanics,
		University of Science and Technology of China, Hefei.
		zimoliao@mail.ustc.edu.cn
*/
class LatticeSolver
{
	/* parameters */
	// MPI
	const int host = 0;
	int size, rank;

	// geometry
	int nx, ny, di, dj, b[9]; // global
	int ni, nj; // local

	// model
	bool is_srt, is_mrt;
	double tau, omega_e, omega_ep, omega_q, omega_nu;

	// file stream
	ifstream fin;
	ofstream fout;


	/* data */
	int step = 0;

	// moment & population
	LatticeMoment lm;
	LatticePopulation lp;

	// fluid-boundaries & entities (particles)
	int nlb, ns, ne;
	int nle;
	vector<LatticeEntity> le;


	/* internal functions */
	void Ignore(int l);

	void InitParameter();
	void InitBoundary();
	void InitEntity();

	void StartOrder();
	void EndOrder();

public:
	/* constructor */
	LatticeSolver();

	/* control functions */
	void Calculate();

	void PrintInfo();

	void WriteUnit();
	void WriteFlow();

=======
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

>>>>>>> new
};


#endif // !LATTICESOLVER_H_