#ifndef LATTICEBOUND_H_
#define LATTICEBOUND_H_

#include<algorithm>
#include<iostream> // TODO DELETE

using std::max;
using std::min;

class LatticeBound
{
	/* parameters */
	int type; // 1:NEBB-0; 2:NEBB-V; 3:NEBB-P;
	int direct;
	double rho, u, v;

	bool exist;
	int imin, imax, jmin, jmax;
	int istart, iend, jstart, jend;

public:


	/* constructor */
	LatticeBound(int rank, int ni, int nj, \
		int type, int direct, \
		int istart, int iend, int jstart, int jend, \
		double rho, double u, double v);


	/* functions */
	bool IsExist();

	void CalculateNebb0(double* f);
	void CalculateNebb0(double* f, double rho_ref);
	void CalculateNebbV(double* f);
	void CalculateNebbP(double* f);

	int GetType();
	int GetIs();
	int GetIe();
	int GetJs();
	int GetJe();
};


#endif // !LATTICEBOUND_H_