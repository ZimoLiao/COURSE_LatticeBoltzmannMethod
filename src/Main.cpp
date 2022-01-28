#include<mpi.h>
<<<<<<< HEAD
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
=======
#include<string>

#include"LatticePopulation.h"
#include"LatticeMoment.h"
#include"ImmersedParticle.h"
#include"LatticeBound.h"

#include"LatticeSolver.h"

using std::string;

constexpr int host = 0;

int main()
{
	int rank, size;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	LatticeSolver ls;

	ls.Calculate();

	MPI_Finalize();
	return 0;
>>>>>>> new
}