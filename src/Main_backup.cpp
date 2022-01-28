#include<mpi.h>
#include<string>

#include"LatticePopulation.h"
#include"LatticeMoment.h"
#include"ImmersedParticle.h"
#include"LatticeBound.h"

using std::string;

constexpr int host = 0;

int main_backup()
{
	int rank, size;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	/* configuration */
	int tg[5] = { 0,1,-1,1,-1 };
	if (rank == 0) {
		tg[3] = 1;
	}
	else if (rank == size - 1) {
		tg[1] = 1;
	}

	int ni = 200, nj = 160;
	int nforce = 10;

	int nstep = 0;

	/* initialization */
	LatticeMoment lm;
	lm.InitGeom(rank*ni, 0, ni, nj, tg);
	lm.InitData();
	//lm.InitDataShear();
	lm.UpdateGhost();

	LatticePopulation lp;
	lp.InitGeom(rank*ni, 0, ni, nj, tg);
	lp.InitParam(0.53);
	lp.InitData(lm);
	lp.UpdateGhost();

	ImmersedParticle part(double(rank*ni), 0.0, double(ni), double(nj), \
		25.0, 75.0, 78.0, 0.0, 0.0, 0.0, 0.0);

	int xs = 0, xe = ni * size - 1, ys = 0, ye = nj - 1;

	LatticeBound lb4(rank, ni, nj, 0, 4, xs + 1, xe - 1, ys, ys, 0., 0., 0.);
	LatticeBound lb2(rank, ni, nj, 0, 2, xs + 1, xe - 1, ye, ye, 0., 0., 0.);

	LatticeBound lb3(rank, ni, nj, 2, 3, xe, xe, ys + 1, ye - 1, 0.999, 0., 0.);
	LatticeBound lb1(rank, ni, nj, 2, 1, xs, xs, ys + 1, ye - 1, 1.001, 0., 0.);

	LatticeBound lb5(rank, ni, nj, 0, 7, xe, xe, ye, ye, 0.999, 0., 0.);
	LatticeBound lb6(rank, ni, nj, 0, 8, xs, xs, ye, ye, 1.001, 0., 0.);
	LatticeBound lb7(rank, ni, nj, 0, 5, xs, xs, ys, ys, 1.001, 0., 0.);
	LatticeBound lb8(rank, ni, nj, 0, 6, xe, xe, ys, ys, 0.999, 0., 0.);

	lp.InitBoundary(lb1);
	lp.InitBoundary(lb2);
	lp.InitBoundary(lb3);
	lp.InitBoundary(lb4);
	lp.InitBoundary(lb5);
	lp.InitBoundary(lb6);
	lp.InitBoundary(lb7);
	lp.InitBoundary(lb8);

	/* calculation */
	lm.WriteAscii(0);
	for (int step = 1; step != nstep + 1; step++) {

		for (int i = 0; i != nforce; i++) {
			part.ForceMoment(lm);
		}

		lp.CollideSrt(lm);
		lp.UpdateGhost();
		lp.Stream();
		lp.Boundary();

		lm.Update(lp);

		if (step % 50 == 0) {
			if (rank == host) {
				cout << "step " << step << endl;
			}
			lm.WriteAscii(step);
		}

	}

	MPI_Finalize();
	return 0;
}