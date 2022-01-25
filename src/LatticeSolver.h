#ifndef LATTICESOLVER_H_
#define LATTICESOLVER_H_

#include<mpi.h>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<string>

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

	

public:
	/* constructor */
	LatticeSolver();

	/* control functions */
	void WriteUnit();
	void WriteFlow();

};


#endif // !LATTICESOLVER_H_