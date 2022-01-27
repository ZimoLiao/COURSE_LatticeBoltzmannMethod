#include<mpi.h>

#include"LatticeBlock.h"
#include"LatticeMoment.h"

int main()
{
	int rank, size;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int tg[5] = { 0,1,0,1,0 };
	if (rank == 0) {
		tg[3] = 0;
	}
	else if (rank == size - 1) {
		tg[1] = 0;
	}

	LatticeMoment lm;
	lm.InitGeom(0, 0, 200, 200, tg);
	lm.UpdateGhost();

	MPI_Finalize();

	return 0;
}