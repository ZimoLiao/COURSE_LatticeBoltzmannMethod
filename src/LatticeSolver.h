#ifndef LATTICESOLVER_H_
#define LATTICESOLVER_H_

#include<string>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<mpi.h>

#include"LatticeMoment.h"
#include"LatticePopulation.h"

using namespace std;

constexpr int host = 0;

class LatticeSolver
{
	/* MPI information */
	int size, rank;

	/* simulation parameters */
	// case name
	string case_name;

	// grid
	int Nx, Ny;
	// domain decomposition
	int Dx, Dy;
	// boundary types
	int B[9] = { 0 };

	// collision operators and corresponding relaxation times
	bool SRT = true, MRT = false;
	double tau, omega_e, omega_ep, omega_q, omega_nu;

	// file stream
	fstream fin;

	/* data */
	LatticeMoment lm;
	LatticePopulation lp;

	/* functions */
	void Ignore(int l);

public:
	LatticeSolver();

};

#endif // !LATTICESOLVER_H_