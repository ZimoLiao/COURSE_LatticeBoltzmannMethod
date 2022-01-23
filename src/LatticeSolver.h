#ifndef LATTICESOLVER_H_
#define LATTICESOLVER_H_

#include<fstream>;

#include"LatticeBound.h"
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
	/* global parameters */
	// MPI
	int size, rank;

	// geometry
	int ni, nj, di, dj, b[9];

	// model
	bool is_srt, is_mrt;
	double tau, omega_e, omega_ep, omega_q, omega_nu;

	// file stream
	fstream fin, fout;

	/* data */
	LatticeMoment lm;
	LatticePopulation lp;

	/* internal functions */


public:
	// constructor - mission initialization
	LatticeSolver();

	/* functions */
	// calculation
	void Calculate();
	
	// output
	void OutputAscii();

};

#endif // !LATTICESOLVER_H_