#include "LatticePopulation.h"

<<<<<<< HEAD
inline int LatticePopulation::Index(int i, int j)
{
	return 9 * (nj*(i + 1) + j + 1); // for outer observers
}

void LatticePopulation::CalculateFeq(double * f, double * m)
{
	double udotc, expr = 1.0 - 1.5 * (m[1] * m[1] + m[2] * m[2]);
	f[0] = m[0] * expr * w[0];
	for (int a = 1; a != 9; a++) {
		udotc = m[1] * cx[a] + m[2] * cy[a];
		f[a] = m[0] * (expr + udotc * (3.0 + 4.5 * udotc)) * w[a];
	}
}
=======
>>>>>>> new

void LatticePopulation::InitData(LatticeMoment & lm)
{
<<<<<<< HEAD
	ni = 1;
	nj = 1;
	sizei = 3;
	sizej = 3;
	sizeij = sizei * sizej;
	size = 9 * sizeij;

	data = new double[size];
	send = new double[sizej * 9];
	recv = new double[sizej * 9];

	double feq[9], m[3] = { 1.0,0.0,0.0 };
	CalculateFeq(feq, m);
	for (int i = 0; i != sizeij; i++) {
		data[i] = feq[i % 9];
=======
	FlushData();
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			CalculateEquilibrium(&data[Index(i, j)], &lm(i, j));
		}
>>>>>>> new
	}
}

void LatticePopulation::InitBoundary(LatticeBound & newlb)
{
<<<<<<< HEAD
	delete[] data;
	delete[] send;
	delete[] recv;
}

void LatticePopulation::InitBoundary(LatticeBound newlb)
{
	lb.push_back(newlb);
	nlb++;
}

double & LatticePopulation::operator()(int i, int j)
{
	return data[Index(i, j)];
=======
	if (newlb.IsExist()) {
		lb.push_back(newlb);
		nlb++;
	}
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
	for (int b = 0; b != nlb; b++) {
		switch (lb[b].GetType())
		{
		case 1: // NEBB-0
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					lb[b].CalculateNebb0(&data[Index(i, j)]);
				}
			}
			break;
		case 2: // NEBB-V
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					lb[b].CalculateNebbV(&data[Index(i, j)]);
				}
			}
			break;
		case 3: // NEBB-P
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					lb[b].CalculateNebbP(&data[Index(i, j)]);
				}
			}
			break;
		}
	}
}

void LatticePopulation::CollideSrt(LatticeMoment & lm)
{
	for (int i = 0; i != size; i++) {
		data[i] *= omegac;
	}
	double feq[9];
	int ind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			CalculateEquilibrium(feq, &lm(i, j));
			ind = Index(i, j);
			for (int a = 0; a != 9; a++) {
				data[ind + a] += feq[a] * omega;
			}
		}
	}
>>>>>>> new
}
