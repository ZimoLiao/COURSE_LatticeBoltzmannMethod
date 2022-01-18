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

	/* MPI connection */
	int rank = 0;
	//	0	1	2	3	4	5	6	7	8
	int rank_conn[9] = { rank,rank,rank,rank,rank,rank,rank,rank,rank }; // default to periodic

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

	void InitConnection(); // TODO: initialize the connection information

	/* operators */
	// index	0:n-1
	double& operator()(int i, int j, int p);

	void operator=(const LatticePopulation& lp);

	/* functions */
	void Stream();
	void CollideSrt(LatticeMoment& lm);
	void CollideMrt(LatticeMoment& lm);
	void UpdateGhost();
};

#endif // !LATTICEPOPULATION_H_