#include<mpi.h>
#include<iostream>

#include"LatticePopulation.h"
#include"LatticeMoment.h"

using namespace std;

int main()
{
	MPI_Init(NULL, NULL);

	LatticeMoment lm;
	lm.Init(0, 0, 10, 10);
	LatticePopulation lp;
	lp.Init(lm);
	lp.Stream();


	MPI_Finalize();
	return 0;
}