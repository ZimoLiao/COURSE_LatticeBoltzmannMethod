#include"Head.h"

#ifdef PARALLEL
#include<mpi.h>
#endif // PARALLEL

#include"LatticeBlock.h"

int main()
{
#ifdef PARALLEL
	MPI_Init(NULL, NULL);
#endif // PARALLEL

	LatticeBlock blk;

	int tg[5] = { 0,0,0,0,0 };
	blk.InitGeom(0, 0, 10, 10, tg);

	blk.UpdateGhost();

	cout << "a" << endl;

#ifdef PARALLEL
	MPI_Finalize();
#endif // PARALLEL

	
	return 0;
}