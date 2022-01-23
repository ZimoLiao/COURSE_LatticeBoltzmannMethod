#include "LatticePopulation.h"

inline int LatticePopulation::IndexF(int i, int j)
{
	return 9 * (sizej * i + j);
}

void LatticePopulation::CalculateFeq(double * f, const double * m)
{
	double udotc, expr = 1.0 - 1.5 * (m[1] * m[1] + m[2] * m[2]);
	f[0] = m[0] * expr * w[0];
	for (int a = 1; a != 9; a++) {
		udotc = m[1] * cx[a] + m[2] * cy[a];
		f[a] = m[0] * (expr + udotc * (3.0 + 4.5 * udotc)) * w[a];
	}
}

void LatticePopulation::CalculateFcollMrt(double * f, const double * moment)
{
	double rho = moment[0], u = moment[1], v = moment[2];
	double m[9] = { 0.0 };
	for (int i = 0; i != 9; i++) {
		for (int j = 0; j != 9; j++) {
			m[i] += M[i][j] * f[j];
		}
	}

	double jx = rho * u, jy = rho * v;
	// ref: https://zhuanlan.zhihu.com/p/457061066
	double meq[9] = {
			rho,
			rho * (-2. / 3. + 3. * u * u + 3. * v * v),
			rho * (1. / 3. - 3. * u * u - 3. * v * v) / 3.,
			jx,
			-jx / 3.,
			jy,
			-jy / 3.,
			rho * (u * u - v * v),
			rho * u * v
	};

	// TODO: need optimization!!!
	for (int i = 0; i != 9; i++) { // mstar
		m[i] *= momegac[i];
		m[i] += momega[i] * meq[i];
	}

	for (int i = 0; i != 9; i++) {
		f[i] = 0.0;
		for (int j = 0; j != 9; j++) {
			f[i] += Minv[i][j] * m[j];
		}
	}
}

void LatticePopulation::PackBuffer(double * buf, int direct)
{
	switch (direct)
	{
	case 1:
		for (int i = 0; i != 9; i++) {
			buf[b1 + i] = data[i5 + i];
			buf[b3 + i] = data[i8 + i];
		}
		for (int i = 0; i != 9 * nj; i++) {
			buf[b2 + i] = data[i1 + i];
		}
		break;
	case 3:
		for (int i = 0; i != 9; i++) {
			buf[b1 + i] = data[i6 + i];
			buf[b3 + i] = data[i7 + i];
		}
		for (int i = 0; i != 9 * nj; i++) {
			buf[b2 + i] = data[i3 + i];
		}
		break;
	}
}

void LatticePopulation::UnpackBuffer(double * buf, int direct)
{
	switch (direct)
	{
	case 1:
		for (int i = 0; i != 9 * sizej; i++) {
			data[o8 + i] = buf[i];
		}
		break;
	case 3:
		for (int i = 0; i != 9 * sizej; i++) {
			data[o7 + i] = buf[i];
		}
		break;
	}
}

void LatticePopulation::FlushBuffer(double * buf)
{
	for (int i = 0; i != 9 * sizej; i++) {
		buf[i] = 0.0;
	}
}

LatticePopulation::LatticePopulation()
{
	ni = 1;
	nj = 1;
	sizei = 3;
	sizej = 3;
	sizeij = sizei * sizej;
	size = 9 * sizeij;

	MPI_Comm_rank(MPI_COMM_WORLD, rank);
	for (int i = 1; i != 9; i++) { rank[i] = rank[0]; }

	data = new double[size];
	send = new double[sizej * 9];
	recv = new double[sizej * 9];

	double feq[9], m[3] = { 1.0,0.0,0.0 };
	CalculateFeq(feq, m);
	for (int i = 0; i != sizeij; i++) {
		data[i] = feq[i % 9];
	}
}

LatticePopulation::~LatticePopulation()
{
	delete[] data;
	delete[] send;
	delete[] recv;
}

void LatticePopulation::Init(LatticeMoment & lm)
{
	ni = lm.ni;
	nj = lm.nj;
	sizei = ni + 2;
	sizej = nj + 2;
	sizeij = sizei * sizej;
	size = 9 * sizeij;

	MPI_Comm_rank(MPI_COMM_WORLD, rank);
	for (int i = 1; i != 9; i++) { rank[i] = rank[0]; }

	data = new double[size];
	send = new double[sizej * 9];
	recv = new double[sizej * 9];

	double feq[9];
	int find, mind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			find = IndexF(i + 1, j + 1);
			mind = lm.IndexM(i, j);

			CalculateFeq(feq, &lm.data[mind]);

			for (int a = 0; a != 9; a++) {
				data[find + a] = feq[a];
			}
		}
	}

	// region index
	i1 = IndexF(ni, 1);
	i3 = IndexF(1, 1);
	i5 = IndexF(ni, nj);
	i6 = IndexF(1, nj);
	i7 = IndexF(1, 1);
	i8 = IndexF(ni, 1);

	o1 = IndexF(ni + 1, 1);
	o3 = IndexF(0, 1);
	o5 = IndexF(ni + 1, nj + 1);
	o6 = IndexF(0, nj + 1);
	o7 = IndexF(0, 0);
	o8 = IndexF(ni + 1, 0);

	b1 = 0;
	b2 = 9;
	b3 = 9 * (sizej - 1);
}

void LatticePopulation::InitParameter(double tau)
{
	this->tau = tau;
	omega = 1.0 / tau;
	omegac = 1.0 - omega;
}

void LatticePopulation::InitParameter(double omega_e, double omega_ep, double omega_q, double omega_nu)
{
	this->omega_e = omega_e;
	this->omega_ep = omega_ep;
	this->omega_q = omega_q;
	this->omega_nu = omega_nu;

	momega[1] = omega_e;
	momega[2] = omega_ep;
	momega[3] = omega_q;
	momega[6] = omega_q;
	momega[7] = omega_nu;
	momega[8] = omega_nu;
	for (int i = 0; i != 9; i++) { momegac[i] = 1.0 - momega[i]; }
}

void LatticePopulation::InitConnection(int di, int b[9])
{
	/* series connection topology */

	// connect to other blocks
	rank[2] = rank[0];
	rank[4] = rank[0];

	rank[1] = (rank[0] + 1) % di;
	rank[5] = (rank[0] + 1) % di;
	rank[8] = (rank[0] + 1) % di;

	if (rank[0] != 0) {
		rank[3] = (rank[0] - 1) % di;
		rank[6] = (rank[0] - 1) % di;
		rank[7] = (rank[0] - 1) % di;
	}
	else {
		rank[3] = di - 1;
		rank[6] = di - 1;
		rank[7] = di - 1;
	}

	// extrapolation boundary if specified
	if (rank[0] == 0) {
		if (b[3] < 0) { rank[3] = b[3]; }
		if (b[6] < 0) { rank[6] = b[6]; }
		if (b[7] < 0) { rank[7] = b[7]; }
	}
	if (rank[0] == di - 1) {
		if (b[1] < 0) { rank[1] = b[1]; }
		if (b[5] < 0) { rank[5] = b[5]; }
		if (b[8] < 0) { rank[8] = b[8]; }
	}
	if (b[2] < 0) {
		rank[2] = b[2];
		if (rank[0] != 0) { rank[6] = b[2]; }
		if (rank[0] != di - 1) { rank[5] = b[2]; }
	}
	if (b[4] < 0) {
		rank[4] = b[4];
		if (rank[0] != 0) { rank[7] = b[4]; }
		if (rank[0] != di - 1) { rank[8] = b[4]; }
	}
}

void LatticePopulation::InitBoundary(LatticeBound & lb)
{
	nlb++;
	this->lb.push_back(lb);
}

void LatticePopulation::Stream()
{
	// 1 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 1] = data[9 * (sizej * (i - 1) + j) + 1];
		}
	}

	// 2 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 2] = data[9 * (sizej * (i)+(j - 1)) + 2];
		}
	}

	// 3 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 3] = data[9 * (sizej * (i + 1) + j) + 3];
		}
	}

	// 4 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 4] = data[9 * (sizej * (i)+(j + 1)) + 4];
		}
	}

	// 5 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 5] = data[9 * (sizej * (i - 1) + (j - 1)) + 5];
		}
	}

	// 6 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 6] = data[9 * (sizej * (i + 1) + (j - 1)) + 6];
		}
	}

	// 7 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 7] = data[9 * (sizej * (i + 1) + (j + 1)) + 7];
		}
	}

	// 8 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 8] = data[9 * (sizej * (i - 1) + (j + 1)) + 8];
		}
	}
}

void LatticePopulation::Boundary()
{
	for (int b = 0; b <= nlb; b++) {
		switch (lb[b].type)
		{
		case 0: // NEBB-0
			for (int i = lb[b].is; i != lb[b].ie; i++) {
				for (int j = lb[b].js; j != lb[b].je; j++) {
					lb[b].CalculateNebb0(&data[IndexF(i, j)]);
				}
			}
			break;
		case 1: // NEBB-V
			for (int i = lb[b].is; i != lb[b].ie; i++) {
				for (int j = lb[b].js; j != lb[b].je; j++) {
					lb[b].CalculateNebbV(&data[IndexF(i, j)]);
				}
			}
			break;
		case 2: // NEBB-P
			for (int i = lb[b].is; i != lb[b].ie; i++) {
				for (int j = lb[b].js; j != lb[b].je; j++) {
					lb[b].CalculateNebbP(&data[IndexF(i, j)]);
				}
			}
			break;
		}
	}
}

void LatticePopulation::CollideSrt(LatticeMoment & lm)
{
	for (int i = 0; i != size; i++) { data[i] *= omegac; }

	int find, mind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			find = IndexF(i + 1, j + 1);
			mind = lm.IndexM(i, j);

			double feq[9];
			CalculateFeq(feq, &lm.data[mind]);

			for (int a = 0; a != 9; a++) {
				data[find + a] += feq[a] * omega;
			}
		}
	}
}

void LatticePopulation::CollideMrt(LatticeMoment & lm)
{
	int find, mind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			find = IndexF(i + 1, j + 1);
			mind = lm.IndexM(i, j);

			CalculateFcollMrt(&data[find], &lm.data[mind]);
		}
	}
}

void LatticePopulation::UpdateGhost()
{
	MPI_Status stat;

	// top and bottom sides
	int i0, ig0, j0, jg0;
	if (rank[2] == rank[0]) { // TODO: Optimization!
		jg0 = sizej - 1;
		j0 = 1;
		for (int i = 1; i <= ni; i++) {
			ig0 = 9 * (sizej * i + jg0);
			i0 = 9 * (sizej * i + j0);
			for (int p = 0; p != 9; p++) {
				data[ig0 + p] = data[i0 + p];
			}
		}
	}
	else if (rank[2] < 0) {
		switch (abs(rank[2]))
		{
		default:
			jg0 = sizej - 1;
			j0 = nj;
			for (int i = 1; i <= ni; i++) {
				ig0 = 9 * (sizej * i + jg0);
				i0 = 9 * (sizej * i + j0);
				for (int p = 0; p != 9; p++) {
					data[ig0 + p] = data[i0 + p];
				}
			}
			break;
		}
	}
	if (rank[4] == rank[0]) { // TODO: Optimization!
		jg0 = 0;
		j0 = nj;
		for (int i = 1; i <= ni; i++) {
			ig0 = 9 * (sizej * i + jg0);
			i0 = 9 * (sizej * i + j0);
			for (int p = 0; p != 9; p++) {
				data[ig0 + p] = data[i0 + p];
			}
		}
	}
	else if (rank[4] < 0) {
		switch (abs(rank[4]))
		{
		default:
			jg0 = 0;
			j0 = 1;
			for (int i = 1; i <= ni; i++) {
				ig0 = 9 * (sizej * i + jg0);
				i0 = 9 * (sizej * i + j0);
				for (int p = 0; p != 9; p++) {
					data[ig0 + p] = data[i0 + p];
				}
			}
			break;
		}
	}

	// message-passing
	if (rank[0] == 0) {
		// exchange data with right block
		if (rank[1] >= 0) {
			PackBuffer(send, 1);
			MPI_Send(send, 9 * sizej, MPI_DOUBLE, rank[1], rank[0], MPI_COMM_WORLD);

			MPI_Recv(recv, 9 * sizej, MPI_DOUBLE, rank[1], rank[1], MPI_COMM_WORLD, &stat);
			UnpackBuffer(recv, 1);

#ifdef _DEBUG
			//cout << "message-passing between " << rank[0] << " and " << rank[1] << endl;
#endif // _DEBUG
		}
		if (rank[3] >= 0) {
			// exchange data with left block
			MPI_Recv(recv, 9 * sizej, MPI_DOUBLE, rank[3], rank[3], MPI_COMM_WORLD, &stat);
			UnpackBuffer(recv, 3);

			PackBuffer(send, 3);
			MPI_Send(send, 9 * sizej, MPI_DOUBLE, rank[3], rank[0], MPI_COMM_WORLD);
		}
	}
	else {
		// exchange data with left block
		if (rank[3] >= 0) {
			MPI_Recv(recv, 9 * sizej, MPI_DOUBLE, rank[3], rank[3], MPI_COMM_WORLD, &stat);
			UnpackBuffer(recv, 3);

			PackBuffer(send, 3);
			MPI_Send(send, 9 * sizej, MPI_DOUBLE, rank[3], rank[0], MPI_COMM_WORLD);
		}
		// exchange data with right block
		if (rank[1] >= 0) {
			PackBuffer(send, 1);
			MPI_Send(send, 9 * sizej, MPI_DOUBLE, rank[1], rank[0], MPI_COMM_WORLD);

			MPI_Recv(recv, 9 * sizej, MPI_DOUBLE, rank[1], rank[1], MPI_COMM_WORLD, &stat);
			UnpackBuffer(recv, 1);

#ifdef _DEBUG
			//cout << "message-passing between " << rank[0] << " and " << rank[1] << endl;
#endif // _DEBUG
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	FlushBuffer(send);
	FlushBuffer(recv);


	// extrapolation
	if (rank[5] < 0) {
		switch (abs(rank[5]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data[IndexF(ni + 1, nj + 1) + i] = data[9 * (sizej * ni + nj) + i];
			}
			break;
		}
	}
	if (rank[6] < 0) {
		switch (abs(rank[6]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data[9 * (sizej - 1) + i] = data[9 * (sizej + nj) + i];
			}
			break;
		}
	}
	if (rank[7] < 0) {
		switch (abs(rank[7]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data[i] = data[9 * (sizej + 1) + i];
			}
			break;
		}
	}
	if (rank[8] < 0) {
		switch (abs(rank[8]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data[9 * (sizej * (sizei - 1)) + i] = data[9 * (sizej * ni + 1) + i];
			}
			break;
		}
	}
}

