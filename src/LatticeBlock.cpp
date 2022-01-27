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
				ng[i] = 3;
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
	if (tg[2] < 0) {

	}
	else if (tg[2] > 0) {

	}


	/* left and right sides */
#ifdef PARALLEL
	
	MPI_Status status;

#endif // PARALLEL
	
}
