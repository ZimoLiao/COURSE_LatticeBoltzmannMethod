
#include"LatticeMoment.h"
#include"LatticePopulation.h"

int main()
{

	LatticePopulation lp;
	double d[3] = { 1.0,0.2,-0.05 };
	LatticeMoment lm(20, 20, d);

	lm.OutputAscii();

	lp.Init(lm);

	for (int i = 0; i != 20; i++) {
		lp.CollideSrt(lm);
		lp.UpdateGhost();
		lp.Stream();

		lm.Update(lp);
	}
	lm.OutputAscii();
}