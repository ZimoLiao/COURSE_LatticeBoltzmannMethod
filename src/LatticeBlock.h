#ifndef LATTICEBLOCK_H_
#define LATTICEBLOCK_H_

#include<mpi.h>
#include<iostream>
#include<algorithm>

using std::cout;
using std::endl;
using std::max;

template<const int NV, const int NG>
class LatticeBlock
{
protected:
	const int nvar = NV; // number of variables
	const int ngconn = NG; // number of ghost layers


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

	// pack/unpack buffer-data for message-passing (in x-direction)
	void PackBuffer(double* buf, int direct);
	void UnpackBuffer(double* buf, int direct);

	// calculation
	void CalculateEquilibrium(double* f, const double* m);
	void CalculateMoment(double* m, const double* f);


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


template<const int NV, const int NG>
inline int LatticeBlock<NV, NG>::Index(int i, int j)
{
#ifdef _DEBUG
	if (i < -ng[3] || i >= ni + ng[1] || j < -ng[4] || j >= nj + ng[2]) {
		cout << "err: out_of_range\n";
		return 0;
	}
#endif // _DEBUG

	return nvar * (sizej*(i + ng[3]) + j + ng[4]);
}

template<const int NV, const int NG>
inline int LatticeBlock<NV, NG>::Index(int i, int j, int v)
{
#ifdef _DEBUG
	if (i < -ng[3] || i >= ni + ng[1] || j < -ng[4] || j >= nj + ng[2] \
		|| v < 0 || v >= nvar) {
		cout << "err: out_of_range\n";
		return 0;
	}
#endif // _DEBUG

	return nvar * (sizej*(i + ng[3]) + j + ng[4]) + v;
}

template<const int NV, const int NG>
inline int LatticeBlock<NV, NG>::IndexBuf(int i, int j) // TODO: 可否略去
{
#ifdef _DEBUG
	if (i < 0 || i >= max(ng[1], ng[3]) || j < -ng[4] || j >= nj + ng[2]) { // TODO: 需要检查
		cout << "err: out_of_range\n";
		return 0;
	}
#endif // _DEBUG

	return nvar * (sizej*i + j + ng[4]);
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::FlushData()
{
	for (int i = 0; i != size; i++) {
		data[i] = i;
	}
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::PackBuffer(double * buf, int direct)
{
	int iin, indin0, indbuf0, inddat0;

	switch (direct)
	{
	case 1:
		iin = ni - ng[1];
		break;
	case 3:
		iin = 0;
		break;
	}

	indin0 = Index(iin, -ng[4]);
	for (int ind = 0; ind != nvar * sizej*ng[direct]; ind++) {
		buf[ind] = data[indin0 + ind];
	}
	if (tg[2] > 0) { // periodic in y-direction
		for (int ig = 0; ig != ng[direct]; ig++) {
			inddat0 = Index(iin + ig, -ng[4]);
			indbuf0 = IndexBuf(ig, nj);

			for (int i = 0; i != nvar; i++) {
				buf[indbuf0 + i] = data[inddat0 + i];
			}
		}
	}
	if (tg[4] > 0) { // periodic in y-direction
		for (int ig = 0; ig != ng[direct]; ig++) {
			inddat0 = Index(iin + ig, nj);
			indbuf0 = IndexBuf(ig, -ng[4]);

			for (int i = 0; i != nvar; i++) {
				buf[indbuf0 + i] = data[inddat0 + i];
			}
		}
	}

}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::UnpackBuffer(double * buf, int direct)
{
	int iout, indout0, indbuf0, inddat0;

	switch (direct)
	{
	case 1:
		iout = ni;
		break;
	case 3:
		iout = -ng[3];
		break;
	}

	indout0 = Index(iout, -ng[4]);
	for (int ind = 0; ind != nvar * sizej*ng[direct]; ind++) {
		data[indout0 + ind] = buf[ind];
	}
	if (tg[2] > 0) { // periodic in y-direction
		for (int ig = 0; ig != ng[direct]; ig++) {
			inddat0 = Index(iout + ig, -ng[4]);
			indbuf0 = IndexBuf(ig, nj);

			for (int i = 0; i != nvar; i++) {
				data[inddat0 + i] = buf[indbuf0 + i];
			}
		}
	}
	if (tg[4] > 0) { // periodic in y-direction
		for (int ig = 0; ig != ng[direct]; ig++) {
			inddat0 = Index(iout + ig, nj);
			indbuf0 = IndexBuf(ig, -ng[4]);

			for (int i = 0; i != nvar; i++) {
				data[inddat0 + i] = buf[indbuf0 + i];
			}
		}
	}
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::CalculateEquilibrium(double * f, const double * m)
{
	double udotc, expr = 1.0 - 1.5 * (m[1] * m[1] + m[2] * m[2]);
	f[0] = m[0] * expr * w[0];
	for (int q = 1; q != 9; q++) {
		udotc = m[1] * cx[q] + m[2] * cy[q];
		f[q] = m[0] * (expr + udotc * (3.0 + 4.5 * udotc)) * w[q];
	}
}

template<int NV, int NG>
inline void LatticeBlock<NV, NG>::CalculateMoment(double * m, const double * f)
{
	m[0] = 0.0;
	m[1] = 0.0;
	m[2] = 0.0;

	for (int i = 0; i != 9; i++) {
		m[0] += f[i];
		m[1] += f[i] * cx[i];
		m[2] += f[i] * cy[i];
	}

	m[1] /= m[0];
	m[2] /= m[0];
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::FlushBuffer(double * buf)
{
	for (int i = 0; i != sizebuf; i++) {
		buf[i] = 0.0;
	}
}

template<const int NV, const int NG>
LatticeBlock<NV, NG>::LatticeBlock()
{
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	data = new double[1];
	send = new double[1];
	recv = new double[1];
}

template<const int NV, const int NG>
LatticeBlock<NV, NG>::~LatticeBlock()
{
	delete[] data;
	delete[] send;
	delete[] recv;
}

template<const int NV, const int NG>
double & LatticeBlock<NV, NG>::operator()(int i, int j)
{
	return data[Index(i, j)];
}

template<const int NV, const int NG>
double & LatticeBlock<NV, NG>::operator()(int i, int j, int v)
{
	return data[Index(i, j, v)];
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::InitGeom(int x0, int y0, int ni, int nj, int tg[5])
{
	this->x0 = x0;
	this->y0 = y0;
	this->ni = ni;
	this->nj = nj;
	for (int i = 0; i != 5; i++) {
		this->tg[i] = tg[i];

		if (tg[i] <= 0) { ng[i] = 1; }
		else if (tg[i] == 1) {
			if (i % 2 == 0) {
				ng[i] = 1;
			}
			else {
				ng[i] = ngconn;
			}
		}
	}

	// grid initialization
	sizei = ni + ng[1] + ng[3];
	sizej = nj + ng[2] + ng[4];
	sizeij = sizei * sizej;
	size = sizeij * this->nvar;
	sizebuf = max(ng[1], ng[3])*sizej*nvar;

	// allocation
	data = new double[size];
	send = new double[sizebuf];
	recv = new double[sizebuf];

	FlushData();
	FlushBuffer(send);
	FlushBuffer(recv);
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::InitParam(double tau)
{
	this->tau = tau;
	omega = 1.0 / tau;
	omegac = 1.0 - omega;
}

template<const int NV, const int NG>
void LatticeBlock<NV, NG>::UpdateGhost()
{
	/* top and bottom sides */
	int ind, indg;

	if (tg[2] < 0) { // TODO: direct extrapolation only 
		for (int i = 0; i != ni; i++) {
			ind = Index(i, nj - 1);
			indg = ind + nvar;
			for (int v = 0; v != nvar; v++) {
				data[indg + v] = data[ind + v];
			}
		}
	}
	else if (tg[2] > 0) { // periodic in y-direction
		for (int i = 0; i != ni; i++) {
			ind = Index(i, 0);
			indg = Index(i, nj);
			for (int v = 0; v != nvar; v++) {
				data[indg + v] = data[ind + v];
			}
		}
	}
	if (tg[4] < 0) { // TODO: direct extrapolation only
		for (int i = 0; i != ni; i++) {
			ind = Index(i, 0);
			indg = ind - nvar;
			for (int v = 0; v != nvar; v++) {
				data[indg + v] = data[ind + v];
			}
		}
	}
	else if (tg[4] > 0) { // periodic in y-direction
		for (int i = 0; i != ni; i++) {
			ind = Index(i, nj - 1);
			indg = Index(i, -1);
			for (int v = 0; v != nvar; v++) {
				data[indg + v] = data[ind + v];
			}
		}
	}

	/* left and right sides */
	// extrapolation boundary
	if (tg[1] < 0) {
		ind = Index(ni - 1, -ng[4]);
		indg = Index(ni, -ng[4]);
		for (int i = 0; i != nvar * sizej; i++) {
			data[indg + i] = data[ind + i];
		}
	}
	if (tg[3] < 0) {
		ind = Index(0, -ng[4]);
		indg = Index(-1, -ng[4]);
		for (int i = 0; i != nvar * sizej; i++) {
			data[indg + i] = data[ind + i];
		}
	}

	// connect to other blocks
	if (mpi_size > 1) { // parallel
		MPI_Status status;

		int rank3 = mpi_rank - 1, rank1 = (mpi_rank + 1) % mpi_size;
		if (rank3 < 0) { rank3 += mpi_size; }

		if (mpi_rank % 2 == 0) {
			if (tg[1] > 0) {
				PackBuffer(send, 1);
				MPI_Send(send, sizebuf, MPI_DOUBLE, rank1, mpi_rank, MPI_COMM_WORLD);

				MPI_Recv(recv, sizebuf, MPI_DOUBLE, rank1, rank1, MPI_COMM_WORLD, &status);
				UnpackBuffer(recv, 1);

#ifdef _DEBUG
				//cout << "message-passing between " << mpi_rank << " and " << rank1 << endl;
#endif // _DEBUG
			}

			if (tg[3] > 0) {
				MPI_Recv(recv, sizebuf, MPI_DOUBLE, rank3, rank3, MPI_COMM_WORLD, &status);
				UnpackBuffer(recv, 3);

				PackBuffer(send, 3);
				MPI_Send(send, sizebuf, MPI_DOUBLE, rank3, mpi_rank, MPI_COMM_WORLD);
			}
		}
		else {
			if (tg[3] > 0) {
				MPI_Recv(recv, sizebuf, MPI_DOUBLE, rank3, rank3, MPI_COMM_WORLD, &status);
				UnpackBuffer(recv, 3);

				PackBuffer(send, 3);
				MPI_Send(send, sizebuf, MPI_DOUBLE, rank3, mpi_rank, MPI_COMM_WORLD);
			}

			if (tg[1] > 0) {
				PackBuffer(send, 1);
				MPI_Send(send, sizebuf, MPI_DOUBLE, rank1, mpi_rank, MPI_COMM_WORLD);

				MPI_Recv(recv, sizebuf, MPI_DOUBLE, rank1, rank1, MPI_COMM_WORLD, &status);
				UnpackBuffer(recv, 1);

#ifdef _DEBUG
				//cout << "message-passing between " << mpi_rank << " and " << rank1 << endl;
#endif // _DEBUG
			}
		}
	}
	else { // serial actually

	}
}


#endif // !LATTICEBLOCK_H_