#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

#include<vector>

#include"LatticeEntity.h"

using std::vector;

class LatticeMoment
{

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