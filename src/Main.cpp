#include<mpi.h>
#include<iostream>

#include"LatticeSolver.h"

using namespace std;

int main()
{
	MPI_Init(NULL,NULL);

	LatticeSolver ls;

	ls.PrintInfo();

	ls.Calculate();
	ls.WriteFlow();
	ls.WriteUnit();

	MPI_Finalize();
}