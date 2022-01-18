#ifndef LATTICEPOPULATION_H_
#define LATTICEPOPULATION_H_

#include"LatticeFunction.h"
#include"LatticeMoment.h"

// Populations for D2Q9
class LatticePopulation
{
protected:
	/* grid parameters */
	int ni_, nj_, sizei_, sizej_, sizeij_, size_;

	LatticeFunction func;

	/* data */
	double* data_;

public:
	friend class LatticeMoment;

	/* constructors & destructor */
	LatticePopulation();
	LatticePopulation(int ni, int nj);
	LatticePopulation(int ni, int nj, double rho, double u, double v);
	LatticePopulation(const LatticePopulation& lp);
	LatticePopulation(LatticeMoment& lm);
	~LatticePopulation();

	/* initialization methods */
	void Init(int ni, int nj);
	void Init(int ni, int nj, double rho, double u, double v);
	void Init(const LatticePopulation& lp);
	void Init(LatticeMoment& lm);

	/* operators */
	// index	0:n-1
	double& operator()(int i, int j, int p);

	void operator=(const LatticePopulation& lp);

	/* functions */

};

#endif // !LATTICEPOPULATION_H_