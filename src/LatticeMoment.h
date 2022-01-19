#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

#include<fstream>
#include<iomanip>
#include<string>
#include<iostream>
#include<math.h>

using namespace std;

#include"LatticeFunction.h"
#include"LatticePopulation.h"

// Moments for D2Q9
//	1	rho (density)
//	2	u
//	3	v
class LatticeMoment
{
protected:
	/* parameters */
	int ni_, nj_, sizeij_, size_;

	LatticeFunction func;

	/* data */
	double* data_;
	double diff = 1.0;

public:
	friend class LatticePopulation;

	/* constructors & destructor */
	LatticeMoment();
	LatticeMoment(int ni, int nj);
	LatticeMoment(int ni, int nj, const double d[3]);
	LatticeMoment(const LatticeMoment& lm);
	~LatticeMoment();

	/* initialization methods */
	void Init(int ni, int nj);
	void Init(int ni, int nj, const double d[3]);
	void Init(const LatticeMoment& lm);

	void InitParameter(double tau, double momega1, double momega2, double momega3, double momega4);

	// TODO: test cases
	void SetVelocityShear(double u);

	/* operators */
	// index	0:n-1
	double& operator()(int i, int j, int m);

	void operator=(const LatticeMoment& lm);
	void operator+=(const LatticeMoment& lm);
	void operator-=(const LatticeMoment& lm);
	void operator*=(const LatticeMoment& lm);
	void operator/=(const LatticeMoment& lm);

	/* functions */
	void Update(LatticePopulation& lp);

	/* I/O */
	void OutputAscii();
	void OutputAscii(string fname);
};

#endif // !LATTICEMOMENT_H_