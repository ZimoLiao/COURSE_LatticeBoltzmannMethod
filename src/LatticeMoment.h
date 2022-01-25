#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

#include<vector>

#include"LatticeEntity.h"

using std::vector;

class LatticeMoment
{
	/* parameters */
	// geometry
	int ni, nj, sizeij, size;

	// velocity set
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };

	/* data */
	// moments array
	double* data;

	// entities (solid particles)
	int nle;
	vector<LatticeEntity> le;

public:

	/* constructor & destructor */
	~LatticeMoment();

	/* initialization */
	void InitEntity(LatticeEntity newle);



};


#endif // !LATTICEMOMENT_H_