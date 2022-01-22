#ifndef LATTICEPARAMETER_H_
#define LATTICEPARAMETER_H_

#include<mpi.h>
#include<iostream>

#include"LatticeMoment.h"

using std::cout;
using std::endl;

class LatticePopulation
{
	/* parameters */
	// geometry
	int ni, nj, sizei, sizej, sizeij, size;

	// connection / periodicity
	int rank[9];

	// velocity set
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };
	const double w[9] = {
		4. / 9.,
		1. / 9., 1. / 9., 1. / 9., 1. / 9.,
		1. / 36., 1. / 36., 1. / 36., 1. / 36.
	};

	// model
	double tau, omega, omegac; // SRT

	// region index
	int i1, i3, i5, i6, i7, i8, o1, o3, o5, o6, o7, o8, b1, b2, b3;

	/* data */
	double* data;
	double* send;
	double* recv;
	// TODO: lattice boundary units

	/* internal functions */
	inline int IndexF(int i, int j);

	void CalculateFeq(double* f, const double* m);

	void PackBuffer(double* buf, int direct);
	void UnpackBuffer(double* buf, int direct);
	void FlushBuffer(double* buf);

public:
	friend class LatticeMoment;

	/* constructors & destructor */
	LatticePopulation();
	~LatticePopulation();

	/* initialization */
	void Init(LatticeMoment& lm);

	void InitParameter(double tau); // for SRT
	void InitConnection(int di, int b[9]);

	/* functions */
	void Stream();
	void Boundary();
	void CollideSrt(LatticeMoment& lm);
	void CollideMrt(LatticeMoment& lm);
	void UpdateGhost();
};

#endif // !LATTICEPARAMETER_H_