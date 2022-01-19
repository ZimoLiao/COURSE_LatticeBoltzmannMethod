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
	func.CalculateFeq(feq, rho, u, v);

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

			func.CalculateFeq(feq, &lm(i, j, 0));

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
	func.CalculateFeq(feq, rho, u, v);

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

			func.CalculateFeq(feq, &lm(i, j, 0));

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

void LatticePopulation::InitParameter(double tau, double momega1, double momega2, double momega3, double momega4)
{
	func.Init(tau, momega1, momega2, momega3, momega4);
}

void LatticePopulation::InitConnection(int rank, int Dx, int Dy, int B[9])
{
	this->rank = rank;

	// initialization the connection information
	int di[9], dj[9];
	dj[0] = rank % Dy;
	di[0] = rank / Dy;
	const int ci[9] = { 0,1,0,-1,0,1,-1,-1,1 };
	const int cj[9] = { 0,0,1,0,-1,1,1,-1,-1 };

	bool w, xc, e, s, yc, n;
	for (int i = 0; i != 9; i++) {

		di[i] = (di[0] + ci[i]);
		dj[i] = (dj[0] + cj[i]);

		w = di[i] == -1;
		e = di[i] == Dx;
		xc = !(w || e);
		s = dj[i] == -1;
		n = dj[i] == Dy;
		yc = !(s || n);

		if (e && yc && B[1]) { rank_conn[i] = B[1]; }
		else if (xc && n && B[2]) { rank_conn[i] = B[2]; }
		else if (w && yc && B[3]) { rank_conn[i] = B[3]; }
		else if (xc && s && B[4]) { rank_conn[i] = B[4]; }
		else if (e && n && B[5]) { rank_conn[i] = B[5]; }
		else if (w && n && B[6]) { rank_conn[i] = B[6]; }
		else if (w && s && B[7]) { rank_conn[i] = B[7]; }
		else if (e && s && B[8]) { rank_conn[i] = B[8]; }
		else {
			di[i] = di[i] % Dx;
			dj[i] = dj[i] % Dy;
			if (di[i] < 0) { di[i] += Dx; }
			if (dj[i] < 0) { dj[i] += Dy; }

			rank_conn[i] = Dy * di[i] + dj[i];
		}
	}

	// TODO: 测试部分 待删去
	/*
	if (rank == 3) {
		cout << rank_conn[6] << '\t' << rank_conn[2] << '\t' << rank_conn[5] << '\n';
		cout << rank_conn[3] << '\t' << rank_conn[0] << '\t' << rank_conn[1] << '\n';
		cout << rank_conn[7] << '\t' << rank_conn[4] << '\t' << rank_conn[8] << '\n';
	}*/
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

	// TODO: 把乘更新都内置到func里！
	double omega = func.somega, omegac = func.somegac;
	for (int i = 0; i != size_; i++) { data_[i] *= omegac; }

	double feq[9];
	int ind;
	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			ind = 9 * (sizej_ * (i + 1) + j + 1);

			func.CalculateFeq(feq, &lm(i, j, 0));

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

	int ind;
	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			ind = 9 * (sizej_ * (i + 1) + j + 1);

			func.CalculateFstar(&data_[ind], &lm(i, j, 0));
		}
	}

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
	else if (rank_conn[1] < 0) { // extrapolation boundary
		switch (abs(rank_conn[1]))
		{
		default:
			ig0 = 9 * (sizej_ * (sizei_ - 1) + 1);
			i0 = 9 * (sizej_ * ni_ + 1);
			for (int i = 0; i != 9 * nj_; i++) {
				data_[ig0 + i] = data_[i0 + i];
			}
			break;
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
	else if (rank_conn[2] < 0) {
		switch (abs(rank_conn[2]))
		{
		default:
			jg0 = sizej_ - 1;
			j0 = nj_;
			for (int i = 1; i <= ni_; i++) {
				ig0 = 9 * (sizej_ * i + jg0);
				i0 = 9 * (sizej_ * i + j0);
				for (int p = 0; p != 9; p++) {
					data_[ig0 + p] = data_[i0 + p];
				}
			}
			break;
		}
	}
	if (rank_conn[3] == rank) {
		ig0 = 9;
		i0 = 9 * (sizej_ * ni_ + 1);
		for (int i = 0; i != 9 * nj_; i++) {
			data_[ig0 + i] = data_[i0 + i];
		}
	}
	else if (rank_conn[3] < 0) {
		switch (abs(rank_conn[3]))
		{
		default:
			ig0 = 9;
			i0 = 9 * (sizej_ + 1);
			for (int i = 0; i != 9 * nj_; i++) {
				data_[ig0 + i] = data_[i0 + i];
			}
			break;
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
	else if (rank_conn[4] < 0) {
		switch (abs(rank_conn[4]))
		{
		default:
			jg0 = 0;
			j0 = 1;
			for (int i = 1; i <= ni_; i++) {
				ig0 = 9 * (sizej_ * i + jg0);
				i0 = 9 * (sizej_ * i + j0);
				for (int p = 0; p != 9; p++) {
					data_[ig0 + p] = data_[i0 + p];
				}
			}
			break;
		}

	}

	// corners
	if (rank_conn[5] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[size_ - 9 + i] = data_[9 * (sizej_ + 1) + i];
		}
	}
	else if (rank_conn[5] < 0) {
		switch (abs(rank_conn[5]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data_[size_ - 9 + i] = data_[9 * (sizej_ * ni_ + nj_) + i];
			}
			break;
		}
	}
	if (rank_conn[6] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[9 * (sizej_ - 1) + i] = data_[9 * (sizej_ * ni_ + 1) + i];
		}
	}
	else if (rank_conn[6] < 0) {
		switch (abs(rank_conn[6]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data_[9 * (sizej_ - 1) + i] = data_[9 * (sizej_ + nj_) + i];
			}
			break;
		}
	}
	if (rank_conn[7] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[i] = data_[9 * (sizej_ * ni_ + nj_) + i];
		}
	}
	else if (rank_conn[7] < 0) {
		switch (abs(rank_conn[7]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data_[i] = data_[9 * (sizej_ + 1) + i];
			}
			break;
		}
	}
	if (rank_conn[8] == rank) {
		for (int i = 0; i != 9; i++) {
			data_[9 * (sizej_ * (sizei_ - 1)) + i] = data_[9 * (sizej_ + nj_) + i];
		}
	}
	else if (rank_conn[8] < 0) {
		switch (abs(rank_conn[8]))
		{
		default:
			for (int i = 0; i != 9; i++) {
				data_[9 * (sizej_ * (sizei_ - 1)) + i] = data_[9 * (sizej_ * ni_ + 1) + i];
			}
			break;
		}
	}
}
