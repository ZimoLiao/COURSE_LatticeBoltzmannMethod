#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

#include<fstream>
#include<string>
#include<iomanip>

#include"LatticePopulation.h"

using std::string;
using std::ofstream;
using std::endl;
using std::to_string;

/* Moments <rho, u, v> */
class LatticeMoment
{
	/* parameters */
	// geometry
	int i0, j0, ni, nj, sizeij, size;

	// velocity set
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };


	/* data */
	double* data;


	/* internal functions */
	inline int IndexM(int i, int j);

	void CalculateMoment(double* m, const double* f);

public:
	double diff = 1.0;
	int step = 0;

	friend class LatticePopulation;

	/* constructors & destructor */
	LatticeMoment();
	LatticeMoment(const LatticeMoment& lm);
	~LatticeMoment();

	/* initialization */
	void Init(int i0, int j0, int ni, int nj);

	/* functions */
	void Update(LatticePopulation& lp);

	void OutputAscii(string fname);
	void OutputAscii();

};

#endif // !LATTICEMOMENT_H_