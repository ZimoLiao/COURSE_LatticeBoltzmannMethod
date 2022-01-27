#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

#include<vector>
#include<fstream>
#include<math.h>

#include"LatticeEntity.h"
#include"LatticePopulation.h"

using std::vector;
using std::ofstream;
using std::max;
using std::min;

class LatticeMoment
{
	/* parameters */
	// grid
	int x0, y0, ni, nj, sizeij, size;

	// velocity set
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };


	/* data */
	// moments array
	double* data;
	double diff = 1.0;

	// entities (solid particles)
	int nle;
	vector<LatticeEntity> le;


	/* internal functions */
	inline int Index(int i, int j);

	double FuncD(double r);

	void Calculate(double* m, const double* f);

public:
	friend class LatticePopulation;

	/* constructor & destructor */
	LatticeMoment();
	~LatticeMoment();


	/* operator */
	double& operator()(int i, int j);


	/* initialization */
	void Init(int x0, int y0, int ni, int nj);

	void InitEntity(LatticeEntity newle);


	/* functions */
	void Update(LatticePopulation& lp);

	void Force();

	void WriteAscii(ofstream& fout, int step);
};


#endif // !LATTICEMOMENT_H_