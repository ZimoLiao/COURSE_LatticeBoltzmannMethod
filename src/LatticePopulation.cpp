#include "LatticePopulation.h"

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

LatticePopulation::LatticePopulation()
{
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
	}
}

LatticePopulation::~LatticePopulation()
{
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
}
