
#include"LatticeMoment.h"
#include"LatticePopulation.h"

int main()
{

	LatticePopulation lp;
	double m0[3] = { 1.0,0.0,0.0 };
	LatticeMoment lm(100, 120, m0);

	lm.SetVelocityShear(0.1);
	lm.OutputAscii("output2/output_0.dat");

	lp.Init(lm);

	int step = 0;
	string fname;
	for (int i = 0; i != 4000; i++) {

		step++;

		//lp.CollideSrt(lm);
		lp.CollideMrt(lm);

		lp.UpdateGhost();
		lp.Stream();

		lm.Update(lp);

		if (step % 50 == 0) {
			fname = "output2/output_" + to_string(step) + ".dat";
			lm.OutputAscii(fname);
		}
	}
}