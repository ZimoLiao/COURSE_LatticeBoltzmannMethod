#include<mpi.h>
#include<iostream>

#include"LatticePopulation.h"
#include"LatticeMoment.h"

using namespace std;

int main()
{
	MPI_Init(NULL, NULL);




	MPI_Finalize();
	return 0;
}