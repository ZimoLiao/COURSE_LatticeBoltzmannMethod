#ifndef LATTICEBOUND_H_
#define LATTICEBOUND_H_

/* Fluid boundary unit
	types:
		0	NEBB-0
		1	NEBB-V
		2	NEBB-P
*/
class LatticeBound
{
public:
	int type;

	int direct, is, ie, js, je;

	double rho, u, v;

	/* constructors */
	LatticeBound(int type, int direct, int is, int ie, int js, int je, \
		double rho, double u, double v);

	/* calculation */
	void CalculateNebb0(double* f);
	void CalculateNebbV(double* f);
	void CalculateNebbP(double* f);
};

#endif // !LATTICEBOUND_H_