#include "LatticePopulation.h"

#include<iostream>
using namespace std;

LatticePopulation::LatticePopulation()
{
	ni_ = 1;
	nj_ = 1;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	// equilibrium population with zero velocity
	for (int i = 0; i != sizeij_; i++) {
		data_[9 * i] = 4. / 9.;
		data_[9 * i + 1] = 1. / 9.;
		data_[9 * i + 2] = 1. / 9.;
		data_[9 * i + 3] = 1. / 9.;
		data_[9 * i + 4] = 1. / 9.;
		data_[9 * i + 5] = 1. / 36.;
		data_[9 * i + 6] = 1. / 36.;
		data_[9 * i + 7] = 1. / 36.;
		data_[9 * i + 8] = 1. / 36.;
	}
}

LatticePopulation::LatticePopulation(int ni, int nj)
{
	ni_ = ni;
	nj_ = nj;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	// equilibrium population with zero velocity
	for (int i = 0; i != sizeij_; i++) {

#ifdef _DEBUG
		if (9 * i < 0 || 9 * i + 8 >= size_) {
			cout << "out_of_range\n";
			return;
		}
#endif // _DEBUG

		data_[9 * i] = 4. / 9.;
		data_[9 * i + 1] = 1. / 9.;
		data_[9 * i + 2] = 1. / 9.;
		data_[9 * i + 3] = 1. / 9.;
		data_[9 * i + 4] = 1. / 9.;
		data_[9 * i + 5] = 1. / 36.;
		data_[9 * i + 6] = 1. / 36.;
		data_[9 * i + 7] = 1. / 36.;
		data_[9 * i + 8] = 1. / 36.;
	}
}

LatticePopulation::LatticePopulation(int ni, int nj, double rho, double u, double v)
{
	ni_ = ni;
	nj_ = nj;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	// calculate equilibrium population
	double feq[9];
	func.CalculateEquilibrium(feq, rho, u, v);

	// equilibrium populaiton with uniform specific velocity
	for (int i = 0; i != size_; i++) {
		data_[i] = feq[i % 9];
	}
}

LatticePopulation::LatticePopulation(const LatticePopulation& lp)
{
	ni_ = lp.ni_;
	nj_ = lp.nj_;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = lp.data_[i]; }
}

LatticePopulation::LatticePopulation(LatticeMoment& lm)
{
	ni_ = lm.ni_;
	nj_ = lm.nj_;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	double feq[9];
	int ind;
	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			ind = 9 * (sizej_ * (i + 1) + j + 1);

			func.CalculateEquilibrium(feq, &lm(i, j, 0));

			data_[ind] = feq[0];
			data_[ind + 1] = feq[1];
			data_[ind + 2] = feq[2];
			data_[ind + 3] = feq[3];
			data_[ind + 4] = feq[4];
			data_[ind + 5] = feq[5];
			data_[ind + 6] = feq[6];
			data_[ind + 7] = feq[7];
			data_[ind + 8] = feq[8];
		}
	}
}

LatticePopulation::~LatticePopulation()
{
	delete[] data_;
}

void LatticePopulation::Init(int ni, int nj)
{
	ni_ = ni;
	nj_ = nj;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	// equilibrium population with zero velocity
	for (int i = 0; i != sizeij_; i++) {

#ifdef _DEBUG
		if (9 * i < 0 || 9 * i + 8 >= size_) {
			cout << "out_of_range\n";
			return;
		}
#endif // _DEBUG

		data_[9 * i] = 4. / 9.;
		data_[9 * i + 1] = 1. / 9.;
		data_[9 * i + 2] = 1. / 9.;
		data_[9 * i + 3] = 1. / 9.;
		data_[9 * i + 4] = 1. / 9.;
		data_[9 * i + 5] = 1. / 36.;
		data_[9 * i + 6] = 1. / 36.;
		data_[9 * i + 7] = 1. / 36.;
		data_[9 * i + 8] = 1. / 36.;
	}
}

void LatticePopulation::Init(int ni, int nj, double rho, double u, double v)
{
	ni_ = ni;
	nj_ = nj;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	// calculate equilibrium population
	double feq[9];
	func.CalculateEquilibrium(feq, rho, u, v);

	// equilibrium populaiton with uniform specific velocity
	for (int i = 0; i != size_; i++) {
		data_[i] = feq[i % 9];
	}
}

void LatticePopulation::Init(const LatticePopulation& lp)
{
	ni_ = lp.ni_;
	nj_ = lp.nj_;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = lp.data_[i]; }
}

void LatticePopulation::Init(LatticeMoment& lm)
{
	ni_ = lm.ni_;
	nj_ = lm.nj_;

	sizei_ = ni_ + 2;
	sizej_ = nj_ + 2;
	sizeij_ = sizei_ * sizej_;
	size_ = sizeij_ * 9;
	data_ = new double[size_];

	double feq[9];
	int ind;
	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			ind = 9 * (sizej_ * (i + 1) + j + 1);

			func.CalculateEquilibrium(feq, &lm(i, j, 0));

			data_[ind] = feq[0];
			data_[ind + 1] = feq[1];
			data_[ind + 2] = feq[2];
			data_[ind + 3] = feq[3];
			data_[ind + 4] = feq[4];
			data_[ind + 5] = feq[5];
			data_[ind + 6] = feq[6];
			data_[ind + 7] = feq[7];
			data_[ind + 8] = feq[8];
		}
	}
}

double& LatticePopulation::operator()(int i, int j, int p)
{

#ifdef _DEBUG

	if (i < -1 || j < -1 || p <= -9 || i > ni_ || j > nj_ || p >= 9) {
		cout << "out_of_range\n";
		return data_[0];
	}

#endif // _DEBUG

	// considering ghost layer
	i++;
	j++;

	if (p >= 0) {
		return data_[9 * (sizej_ * i + j) + p];
	}
	else {
		switch (p)
		{
		case -1:
			return data_[9 * (sizej_ * i + j) + 3];
		case -2:
			return data_[9 * (sizej_ * i + j) + 4];
		case -3:
			return data_[9 * (sizej_ * i + j) + 1];
		case -4:
			return data_[9 * (sizej_ * i + j) + 2];
		case -5:
			return data_[9 * (sizej_ * i + j) + 7];
		case -6:
			return data_[9 * (sizej_ * i + j) + 8];
		case -7:
			return data_[9 * (sizej_ * i + j) + 5];
		case -8:
			return data_[9 * (sizej_ * i + j) + 6];
		}
		return data_[0];
	}
}

void LatticePopulation::operator=(const LatticePopulation& lp)
{

#ifdef _DEBUG

	if (ni_ != lp.ni_ || size_ != lp.size_ || sizeij_ != lp.sizeij_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	for (int i = 0; i != size_; i++) { data_[i] = lp.data_[i]; }
}

void LatticePopulation::Stream()
{
	// 1 direction
	for (int i = ni_; i > 0; i--) {
		for (int j = nj_; j > 0; j--) {
			data_[9 * (sizej_ * i + j) + 1] = data_[9 * (sizej_ * (i - 1) + j) + 1];
		}
	}

	// 2 direction
	for (int i = ni_; i > 0; i--) {
		for (int j = nj_; j > 0; j--) {
			data_[9 * (sizej_ * i + j) + 2] = data_[9 * (sizej_ * (i)+(j - 1)) + 2];
		}
	}

	// 3 direction
	for (int i = 1; i <= ni_; i++) {
		for (int j = 1; j <= nj_; j++) {
			data_[9 * (sizej_ * i + j) + 3] = data_[9 * (sizej_ * (i + 1) + j) + 3];
		}
	}

	// 4 direction
	for (int i = 1; i <= ni_; i++) {
		for (int j = 1; j <= nj_; j++) {
			data_[9 * (sizej_ * i + j) + 4] = data_[9 * (sizej_ * (i)+(j + 1)) + 4];
		}
	}

	// 5 direction
	for (int i = ni_; i > 0; i--) {
		for (int j = nj_; j > 0; j--) {
			data_[9 * (sizej_ * i + j) + 5] = data_[9 * (sizej_ * (i - 1) + (j - 1)) + 5];
		}
	}

	// 6 direction
	for (int i = 1; i <= ni_; i++) {
		for (int j = 1; j <= nj_; j++) {
			data_[9 * (sizej_ * i + j) + 6] = data_[9 * (sizej_ * (i + 1) + (j - 1)) + 6];
		}
	}

	// 7 direction
	for (int i = 1; i <= ni_; i++) {
		for (int j = 1; j <= nj_; j++) {
			data_[9 * (sizej_ * i + j) + 7] = data_[9 * (sizej_ * (i + 1) + (j + 1)) + 7];
		}
	}

	// 8 direction
	for (int i = ni_; i > 0; i--) {
		for (int j = nj_; j > 0; j--) {
			data_[9 * (sizej_ * i + j) + 8] = data_[9 * (sizej_ * (i - 1) + (j + 1)) + 8];
		}
	}
}

void LatticePopulation::CollideSrt(LatticeMoment& lm)
{

#ifdef _DEBUG

	if (ni_ != lm.ni_ || nj_ != lm.nj_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	double omega = func.omega, omegac = func.omegac;
	for (int i = 0; i != size_; i++) { data_[i] *= omegac; }

	double feq[9];
	int ind;
	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			ind = 9 * (sizej_ * (i + 1) + j + 1);

			func.CalculateEquilibrium(feq, &lm(i, j, 0));

			data_[ind] += feq[0] * omega;
			data_[ind + 1] += feq[1] * omega;
			data_[ind + 2] += feq[2] * omega;
			data_[ind + 3] += feq[3] * omega;
			data_[ind + 4] += feq[4] * omega;
			data_[ind + 5] += feq[5] * omega;
			data_[ind + 6] += feq[6] * omega;
			data_[ind + 7] += feq[7] * omega;
			data_[ind + 8] += feq[8] * omega;
		}
	}

}

void LatticePopulation::CollideMrt(LatticeMoment& lm)
{

#ifdef _DEBUG

	if (ni_ != lm.ni_ || nj_ != lm.nj_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

}

void LatticePopulation::UpdateGhost()
{
	// TODO: consider the sequence of message-passing

	// 9*(sizej_*(i)+j)


	/* periodic boundaries */
	// sides
	int i0, ig0, j0, jg0;
	if (rank_conn[1] == rank) {
		ig0 = 9 * (sizej_ * (sizei_ - 1) + 1);
		i0 = 9 * (sizej_ + 1);
		for (int i = 0; i != 9 * nj_; i++) {
			data_[ig0 + i] = data_[i0 + i];
		}
	}
	if (rank_conn[2] == rank) { // TODO: Optimization!
		jg0 = sizej_ - 1;
		j0 = 1;
		for (int i = 1; i <= ni_; i++) {
			ig0 = 9 * (sizej_ * i + jg0);
			i0 = 9 * (sizej_ * i + j0);
			for (int p = 0; p != 9; p++) {
				data_[ig0 + p] = data_[i0 + p];
			}
		}
	}
	if (rank_conn[3] == rank) {
		ig0 = 9;
		i0 = 9 * (sizej_ * ni_ + 1);
		for (int i = 0; i != 9 * nj_; i++) {
			data_[ig0 + i] = data_[i0 + i];
		}
	}
	if (rank_conn[4] == rank) { // TODO: Optimization!
		jg0 = 0;
		j0 = nj_;
		for (int i = 1; i <= ni_; i++) {
			ig0 = 9 * (sizej_ * i + jg0);
			i0 = 9 * (sizej_ * i + j0);
			for (int p = 0; p != 9; p++) {
				data_[ig0 + p] = data_[i0 + p];
			}
		}
	}

	// corners
	if (rank_conn[5] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[size_ - 9 + i] = data_[9 * (sizej_ + 1) + i];
		}
	}
	if (rank_conn[6] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[9 * (sizej_ - 1) + i] = data_[9 * (sizej_ * ni_ + 1) + i];
		}
	}
	if (rank_conn[7] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[i] = data_[9 * (sizej_ * ni_ + nj_) + i];
		}
	}
	if (rank_conn[8] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[9 * (sizej_ * (sizei_ - 1)) + i] = data_[9 * (sizej_ + nj_) + i];
		}
	}
}
