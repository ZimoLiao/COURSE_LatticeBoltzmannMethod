#ifndef LATTICEPOPULATION_H_
#define LATTICEPOPULATION_H_

#include<vector>

#include"LatticeBlock.h"
#include"LatticeMoment.h"
#include"LatticeBound.h"

using std::vector;

class LatticePopulation :
	public LatticeBlock<9, 1>
{

	/* data */
	// fluid-boundaries
	int nlb = 0;
	vector<LatticeBound> lb;


public:
	friend class LatticeMoment;


	/* initialization */
	void InitData(LatticeMoment& lm);

	void InitBoundary(LatticeBound& newlb);



	/* functions */
	void Stream();
	void Boundary();
	void CollideSrt(LatticeMoment& lm);

};

#endif // !LATTICEPOPULATION_H_