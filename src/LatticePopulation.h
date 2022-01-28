#ifndef LATTICEPOPULATION_H_
#define LATTICEPOPULATION_H_

<<<<<<< HEAD
#include<mpi.h>
#include<vector>
#include<iostream>

#include"LatticeBound.h"
=======
#include<vector>

#include"LatticeBlock.h"
>>>>>>> new
#include"LatticeMoment.h"
#include"LatticeBound.h"

using std::vector;
<<<<<<< HEAD
using std::cin;
using std::cout;
using std::endl;

class LatticePopulation
{
	/* parameters */
	// grid
	int ni, nj, sizei, sizej, sizeij, size;

	// region index
	int i1, i3, i5, i6, i7, i8, o1, o3, o5, o6, o7, o8, b1, b2, b3;

	// velocity set
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };
	const double w[9] = {
		4. / 9.,
		1. / 9., 1. / 9., 1. / 9., 1. / 9.,
		1. / 36., 1. / 36., 1. / 36., 1. / 36.
	};

	// model (SRT)
	double tau, omega, omegac;



	/* data */
	// populations array
	double* data;

	// buffer for ghost layer passing
	double* send;
	double* recv;

	// fluid-boundaries
	int nlb;
	vector<LatticeBound> lb;


	/* internal functions */
	inline int Index(int i, int j);

	void CalculateFeq(double* f, double* m);
=======

class LatticePopulation :
	public LatticeBlock<9, 1>
{

	/* data */
	// fluid-boundaries
	int nlb = 0;
	vector<LatticeBound> lb;

>>>>>>> new

public:
	friend class LatticeMoment;

<<<<<<< HEAD
	/* constructor & destructor */
	LatticePopulation();
	LatticePopulation(LatticeMoment& lm);
	~LatticePopulation();


	/* initialization */
	void Init();

	void InitParameter(double tau); // for SRT

	void InitBoundary(LatticeBound newlb);


	/* operator */
	double& operator()(int i, int j);


	/* functions */
=======

	/* initialization */
	void InitData(LatticeMoment& lm);

	void InitBoundary(LatticeBound& newlb);



	/* functions */
	void Stream();
	void Boundary();
	void CollideSrt(LatticeMoment& lm);
>>>>>>> new

};


#endif // !LATTICEPOPULATION_H_