#include<mpi.h>
#include<iostream>

#include"LatticeSolver.h"

using namespace std;

int main()
{
	MPI_Init(NULL,NULL);

	LatticeSolver ls;

	MPI_Finalize();
}