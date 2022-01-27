#include "LatticeBlock.h"

inline int LatticeBlock::Index(int i, int j)
{
#ifdef _DEBUG
	if (i < -ng[3] || i >= sizei || j < -ng[4] || j >= sizej) {
		cout << "err: out_or_range\n";
		return 0;
	}
#endif // _DEBUG

	return nvar * (sizej*(i + ng[3]) + j + ng[4]);
}

inline int LatticeBlock::Index(int i, int j, int v)
{
#ifdef _DEBUG
	if (i < -ng[3] || i >= sizei || j < -ng[4] || j >= sizej \
		|| v < 0 || v >= nvar) {
		cout << "err: out_or_range\n";
		return 0;
	}
#endif // _DEBUG

	return nvar * (sizej*(i + ng[3]) + j + ng[4]) + v;
}

inline int LatticeBlock::IndexBuf(int i, int j) // TODO: 可否略去
{
#ifdef _DEBUG
	if (i < 0 || i >= max(ng[1], ng[3]) || j < -ng[4] || j >= sizej) { // TODO: 需要检查
		cout << "err: out_or_range\n";
		return 0;
	}
#endif // _DEBUG

	return nvar * (sizej*i + j + ng[4]);
}

void LatticeBlock::FlushData()
{
	for (int i = 0; i != size; i++) {
		data[i] = i;
	}
}

void LatticeBlock::PackBuffer(double * buf, int direct)
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

void LatticeBlock::UnpackBuffer(double * buf, int direct)
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

void LatticeBlock::FlushBuffer(double * buf)
{
	for (int i = 0; i != sizebuf; i++) {
		buf[i] = 0.0;
	}
}

LatticeBlock::LatticeBlock()
{
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	data = new double[1];
	send = new double[1];
	recv = new double[1];
}

LatticeBlock::~LatticeBlock()
{
	delete[] data;
	delete[] send;
	delete[] recv;
}

double & LatticeBlock::operator()(int i, int j)
{
	return data[Index(i, j)];
}

double & LatticeBlock::operator()(int i, int j, int v)
{
	return data[Index(i, j, v)];
}

void LatticeBlock::InitGeom(int x0, int y0, int ni, int nj, int tg[5])
{
	this->x0 = x0;
	this->y0 = y0;
	this->ni = ni;
	this->nj = nj;
	for (int i = 0; i != 5; i++) {
		this->tg[i] = tg[i];

		if (tg[i] == 0) { ng[i] = 0; }
		else if (tg[i] < 0) { ng[i] = 1; }
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
	size = sizeij * nvar;
	sizebuf = max(ng[1], ng[3])*sizej*nvar;

	// allocation
	data = new double[size];
	send = new double[sizebuf];
	recv = new double[sizebuf];

	FlushData();
	FlushBuffer(send);
	FlushBuffer(recv);
}

void LatticeBlock::InitParam(double tau)
{
	this->tau = tau;
	omega = 1.0 / tau;
	omegac = 1.0 - omega;
}

void LatticeBlock::UpdateGhost()
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
				cout << "message-passing between " << mpi_rank << " and " << rank1 << endl;
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
				cout << "message-passing between " << mpi_rank << " and " << rank1 << endl;
#endif // _DEBUG
			}
		}
	}
	else { // serial actually

	}
}
