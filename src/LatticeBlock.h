#ifndef LATTICEBLOCK_H_
#define LATTICEBLOCK_H_

#include<mpi.h>
#include<iostream>
#include<algorithm>

using std::cout;
using std::endl;
using std::max;

class LatticeBlock
{
protected:
	/* MPI */
	int mpi_rank, mpi_size;


	/* physics parameters */
	// velocity set (D2Q9)
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };
	const double w[9] = { 4. / 9.,1. / 9., 1. / 9., 1. / 9., 1. / 9.,1. / 36., 1. / 36., 1. / 36., 1. / 36. };

	// collision operator (SRT)
	double tau, omega, omegac;


	/* geometry parameters */
	int x0, y0; // origin position in global coordinate

	const int nvar = 1; // number of variables
	const int ngconn = 1; // number of ghost layers

	int ni, nj, sizei, sizej, sizeij, size, sizebuf; // size of data

	int ng[5] = { 0 }; // size of ghost layers
	int tg[5] = { 0 }; // type of ghost layers


	/* data */
	double* data;
	double* send;
	double* recv;


	/* internal functions */
	// index
	inline int Index(int i, int j);
	inline int Index(int i, int j, int v);
	inline int IndexBuf(int i, int j);

	void FlushData();
	void FlushBuffer(double* buf);

	// pack buffer-data for MPI sending
	void PackBuffer(double* buf, int direct);
	// unpack buffer-data received through MPI
	void UnpackBuffer(double* buf, int direct);


public:
	/* constructors & destructor */
	LatticeBlock();
	~LatticeBlock();


	/* operators */
	double& operator()(int i, int j);
	double& operator()(int i, int j, int v);


	/* initialization */
	void InitGeom(int x0, int y0, int ni, int nj, int tg[5]);
	void InitParam(double tau); // for SRT

	void UpdateGhost();

};


#endif // !LATTICEBLOCK_H_