
#include"LatticeMoment.h"
#include"LatticePopulation.h"

int main()
{

	LatticePopulation lp;
	LatticeMoment lm(100, 100);

	lm.SetVelocityShear(0.1);
	lm.OutputAscii("output_0.dat");

	lp.Init(lm);

	int step = 0;
	string fname;
	for (int i = 0; i != 1000; i++) {

		step++;

		lp.CollideSrt(lm);
		lp.UpdateGhost();
		lp.Stream();

		lm.Update(lp);

		if (step % 50 == 0) {
			fname = "output_" + to_string(step) + ".dat";
			lm.OutputAscii(fname);
		}
	}
}