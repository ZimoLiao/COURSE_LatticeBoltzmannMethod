#ifndef LATTICEBOUND_H_
#define LATTICEBOUND_H_

#include<algorithm>

class LatticeBound
{
	/* index */
	int index_global;

	/* parameters */
	int direct;
	double rho, u, v; // for NEBB boundaries

public:
	/* parameters */
	int type; // 1:NEBB-0; 2:NEBB-V; 3:NEBB-P
	int is, ie, js, je;


	/* constructors */
	LatticeBound(int index_global, int type, int direct, \
		int is, int ie, int js, int je, double rho, double u, double v);

	void UpdateLocal(int rank, int ni);

	bool IsExist(int rank, int ni);

	/* calculation */
	void CalculateNebb0(double* f);
	void CalculateNebbV(double* f);
	void CalculateNebbP(double* f);
};


#endif // !LATTICEBOUND_H_