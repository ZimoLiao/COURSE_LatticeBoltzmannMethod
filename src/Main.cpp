#include<mpi.h>
#include<iostream>

#include"LatticeSolver.h"

using namespace std;

int main()
{
	MPI_Init(NULL, NULL);

	LatticeSolver ls;

	MPI_Barrier(MPI_COMM_WORLD);

	ls.OutputAscii();
	int step = 0;
	for (int i = 0; i != 2000; i++) {
		ls.Calculate();
		step++;

		if (step % 100 == 0) {
			ls.OutputAscii();
		}
	}

	MPI_Finalize();
	return 0;
}