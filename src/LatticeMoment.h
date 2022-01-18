#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

// Moments for D2Q9
//	1	rho (density)
//	2	u
//	3	v
class LatticeMoment
{
protected:
	/* parameters */
	int ni_, nj_, sizeij_, size_;

	/* data */
	double* data_;

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

	/* operators */
	// index	0:n-1
	double& operator()(int i, int j, int m);

	void operator=(const LatticeMoment& lm);
	void operator+=(const LatticeMoment& lm);
	void operator-=(const LatticeMoment& lm);
	void operator*=(const LatticeMoment& lm);
	void operator/=(const LatticeMoment& lm);
};

#endif // !LATTICEMOMENT_H_