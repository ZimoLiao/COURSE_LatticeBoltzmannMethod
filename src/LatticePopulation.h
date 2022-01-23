#ifndef LATTICEPARAMETER_H_
#define LATTICEPARAMETER_H_

#include<mpi.h>
#include<iostream>
#include<vector>

#include"LatticeMoment.h"
#include"LatticeBound.h"

using std::cout;
using std::endl;
using std::vector;

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
	double omega_e, omega_ep, omega_q, omega_nu;
	double momega[9] = { 0.0 };
	double momegac[9] = { 1.0 };
	const double M[9][9] = {
		{	1,	1,	1,	1,	1,	1,	1,	1,	1	},
		{	-4,	-1,	-1,	-1,	-1,	2,	2,	2,	2	},
		{	4,	-2,	-2,	-2,	-2,	1,	1,	1,	1	},
		{	0,	1,	0,	-1,	0,	1,	-1,	-1,	1	},
		{	0,	-2,	0,	2,	0,	1,	-1,	-1,	1	},
		{	0,	0,	1,	0,	-1,	1,	1,	-1,	-1	},
		{	0,	0,	-2,	0,	2,	1,	1,	-1,	-1	},
		{	0,	1,	-1,	1,	-1,	0,	0,	0,	0	},
		{	0,	0,	0,	0,	0,	1,	-1,	1,	-1	},
	};
	const double Minv[9][9] = {
		{	1. / 9,	-1. / 9,	1. / 9,		0,			0,			0,			0,			0,			0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	1. / 6,		-1. / 6,	0,			0,			1. / 4,		0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	0,			0,			1. / 6,		-1. / 6,	-1. / 4,	0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	-1. / 6,	1. / 6,		0,			0,			1. / 4,		0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	0,			0,			-1. / 6,	1. / 6,		-1. / 4,	0		},
		{	1. / 9,	1. / 18,	1. / 36,	1. / 6,		1. / 12,	1. / 6,		1. / 12,	0,			1. / 4	},
		{	1. / 9,	1. / 18,	1. / 36,	-1. / 6,	-1. / 12,	1. / 6,		1. / 12,	0,			-1. / 4	},
		{	1. / 9,	1. / 18,	1. / 36,	-1. / 6,	-1. / 12,	-1. / 6,	-1. / 12,	0,			1. / 4	},
		{	1. / 9,	1. / 18,	1. / 36,	1. / 6,		1. / 12,	-1. / 6,	-1. / 12,	0,			-1. / 4	},
	};

	// region index
	int i1, i3, i5, i6, i7, i8, o1, o3, o5, o6, o7, o8, b1, b2, b3;

	/* data */
	double* data;
	double* send;
	double* recv;
	// TODO: lattice boundary units
	vector<LatticeBound> lb;
	int nlb = -1;

	/* internal functions */
	inline int IndexF(int i, int j);

	void CalculateFeq(double* f, const double* m);
	void CalculateFcollMrt(double* f, const double* moment);

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
	void InitParameter(double omega_e, double omega_ep, double omega_q, double omega_nu); // for MRT

	void InitConnection(int di, int b[9]);
	void InitBoundary(LatticeBound& lb);

	/* functions */
	void Stream();
	void Boundary();
	void CollideSrt(LatticeMoment& lm);
	void CollideMrt(LatticeMoment& lm);
	void UpdateGhost();
};

#endif // !LATTICEPARAMETER_H_