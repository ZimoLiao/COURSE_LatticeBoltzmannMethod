
#include"LatticeMoment.h"
#include"LatticePopulation.h"

int main()
{

	LatticePopulation lp;
	double d[3] = { 1.0,0.2,-0.1 };
	LatticeMoment lm(10, 10, d);

	lm.OutputAscii();

	lp.Init(lm);
	lm.Update(lp);
	lm.OutputAscii();
}