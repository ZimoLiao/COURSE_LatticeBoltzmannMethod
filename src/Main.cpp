#include<mpi.h>
#include<iostream>

#include"LatticeMoment.h"
#include"LatticePopulation.h"
#include"LatticeSolver.h"

using namespace std;

int main()
{
	constexpr int host = 0;
	int size, rank;
	MPI_Status status;

	// initialization
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	LatticeSolver ls;

	ls.Update();

	// sequence computing test
	/*
	LatticePopulation lp;
	double m0[3] = { 1.0,0.0,0.0 };
	LatticeMoment lm(20, 20, m0);

	lm.SetVelocityShear(0.1);
	//lm.OutputAscii("output2/output_0.dat");

	lp.Init(lm);

	int step = 0;
	string fname;
	for (int i = 0; i != 0; i++) {

		step++;

		//lp.CollideSrt(lm);
		lp.CollideMrt(lm);

		lp.UpdateGhost();
		lp.Stream();

		lm.Update(lp);

		if (step % 10 == 0) {
			cout << "step " << step << endl;

			//fname = "output2/output_" + to_string(step) + ".dat";
			//lm.OutputAscii(fname);
		}
	}
	*/


	MPI_Finalize();
	return 0;
}