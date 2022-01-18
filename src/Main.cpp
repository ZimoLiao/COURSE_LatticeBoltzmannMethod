
#include"LatticeMoment.h"
#include"LatticePopulation.h"

int main()
{

	LatticePopulation lp;
	
	LatticeMoment lm(10, 10);

	lm.Init(2, 10);

	lp.Init(lm);

	LatticePopulation lp2;

	lp2.Init(lp);

}