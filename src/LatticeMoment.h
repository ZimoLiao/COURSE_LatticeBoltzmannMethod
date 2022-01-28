#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

<<<<<<< HEAD
#include<vector>
#include<fstream>
#include<math.h>

#include"LatticeEntity.h"
=======
#include<string>
#include<fstream>
#include<iomanip>

#include"LatticeBlock.h"
>>>>>>> new
#include"LatticePopulation.h"
#include"ImmersedParticle.h"

<<<<<<< HEAD
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
=======
using std::string;

constexpr double pi = 3.1415926535897932;

class LatticeMoment :
	public LatticeBlock<3, 3>
{


>>>>>>> new

	// entities (solid particles)
	int nle;
	vector<LatticeEntity> le;


	/* internal functions */
	inline int Index(int i, int j);

	double FuncD(double r);

	void Calculate(double* m, const double* f);

public:
	friend class LatticePopulation;
	friend class ImmersedParticle;

<<<<<<< HEAD
	/* constructor & destructor */
	LatticeMoment();
	~LatticeMoment();


	/* operator */
	double& operator()(int i, int j);


	/* initialization */
	void Init(int x0, int y0, int ni, int nj);

	void InitEntity(LatticeEntity newle);

=======

	/* monitors */
	double diff = 1.0;


	/* initialization */
	void InitData();

	// TODO: just for test
	void InitDataShear();
	// TODO: need delete

>>>>>>> new

	/* functions */
	void Update(LatticePopulation& lp);

<<<<<<< HEAD
	void Force();

	void WriteAscii(ofstream& fout, int step);
=======
	void WriteAscii(int step);

>>>>>>> new
};


#endif // !LATTICEMOMENT_H_