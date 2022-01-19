#include "LatticeMoment.h"

LatticeMoment::LatticeMoment()
{
	ni_ = 1;
	nj_ = 1;

	sizeij_ = 1;
	size_ = 1;
	data_ = new double[1];
	data_[0] = 0.0;
}

LatticeMoment::LatticeMoment(int ni, int nj)
{
	ni_ = ni;
	nj_ = nj;

	sizeij_ = ni_ * nj_;
	size_ = sizeij_ * 3;
	data_ = new double[size_];
	for (int i = 0; i != sizeij_; i++) {

#ifdef _DEBUG
		if (3 * i < 0 || 3 * i + 2 >= size_) {
			cout << "out_of_range\n";
			return;
		}
#endif // _DEBUG

		data_[3 * i] = 1.0;
		data_[3 * i + 1] = 0.0;
		data_[3 * i + 2] = 0.0;
	}
}

LatticeMoment::LatticeMoment(int ni, int nj, const double d[3])
{
	ni_ = ni;
	nj_ = nj;

	sizeij_ = ni_ * nj_;
	size_ = sizeij_ * 3;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = d[i % 3]; }
}

LatticeMoment::LatticeMoment(const LatticeMoment& lm)
{
	ni_ = lm.ni_;
	nj_ = lm.nj_;

	sizeij_ = lm.sizeij_;
	size_ = lm.size_;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = lm.data_[i]; }
}

LatticeMoment::~LatticeMoment()
{
	delete[] data_;
}

void LatticeMoment::Init(int ni, int nj)
{
	ni_ = ni;
	nj_ = nj;

	sizeij_ = ni_ * nj_;
	size_ = sizeij_ * 3;
	data_ = new double[size_];
	for (int i = 0; i != sizeij_; i++) {

#ifdef _DEBUG
		if (3 * i < 0 || 3 * i + 2 >= size_) {
			cout << "out_of_range\n";
			return;
		}
#endif // _DEBUG

		data_[3 * i] = 1.0;
		data_[3 * i + 1] = 0.0;
		data_[3 * i + 2] = 0.0;
	}
}

void LatticeMoment::Init(int ni, int nj, const double d[3])
{
	ni_ = ni;
	nj_ = nj;

	sizeij_ = ni_ * nj_;
	size_ = sizeij_ * 3;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = d[i % 3]; }
}

void LatticeMoment::Init(const LatticeMoment& lm)
{
	ni_ = lm.ni_;
	nj_ = lm.nj_;

	sizeij_ = lm.sizeij_;
	size_ = lm.size_;
	data_ = new double[size_];
	for (int i = 0; i != size_; i++) { data_[i] = lm.data_[i]; }
}

void LatticeMoment::SetVelocityShear(double u)
{
	double nmode = 4.0;

	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			if (j > nj_ / 2 - 1) {
				data_[3 * (nj_ * i + j) + 1] = u;
			}
			else {
				data_[3 * (nj_ * i + j) + 1] = -u;
			}

			// disturbance
			if (j > nj_ * (nmode - 0.5) / 2 / nmode - 1 && j < nj_ * (nmode + 0.5) / 2 / nmode - 1) {
				data_[3 * (nj_ * i + j) + 2] = 0.01 * sin(6.283185307 * (nmode * i / ni_))\
					* cos(3.141592654 * (j - nj_ / 2. + 1.) / (nj_ / 2 / nmode));
			}
		}
	}
}

double& LatticeMoment::operator()(int i, int j, int m)
{

#ifdef _DEBUG

	if (i < 0 || j < 0 || m < 0 || i >= ni_ || j >= nj_ || m >= 3) {
		cout << "out_of_range\n";
		return data_[0];
	}

#endif // _DEBUG

	return data_[3 * (nj_ * i + j) + m];
}

void LatticeMoment::operator=(const LatticeMoment& lm)
{

#ifdef _DEBUG

	if (size_ != lm.size_ || sizeij_ != lm.sizeij_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	for (int i = 0; i != size_; i++) { data_[i] = lm.data_[i]; }
}

void LatticeMoment::operator+=(const LatticeMoment& lm)
{

#ifdef _DEBUG

	if (size_ != lm.size_ || sizeij_ != lm.sizeij_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	for (int i = 0; i != size_; i++) { data_[i] += lm.data_[i]; }
}

void LatticeMoment::operator-=(const LatticeMoment& lm)
{

#ifdef _DEBUG

	if (size_ != lm.size_ || sizeij_ != lm.sizeij_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	for (int i = 0; i != size_; i++) { data_[i] -= lm.data_[i]; }
}

void LatticeMoment::operator*=(const LatticeMoment& lm)
{

#ifdef _DEBUG

	if (size_ != lm.size_ || sizeij_ != lm.sizeij_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	for (int i = 0; i != size_; i++) { data_[i] *= lm.data_[i]; }
}

void LatticeMoment::operator/=(const LatticeMoment& lm)
{

#ifdef _DEBUG

	if (size_ != lm.size_ || sizeij_ != lm.sizeij_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	for (int i = 0; i != size_; i++) { data_[i] /= lm.data_[i]; }
}

void LatticeMoment::Update(LatticePopulation& lp)
{

#ifdef _DEBUG

	if (ni_ != lp.ni_ || nj_ != lp.nj_) {
		cout << "size_mismatch\n";
		return;
	}

#endif // _DEBUG

	diff = 0.0;

	int pind, mind;
	double u, v;
	for (int i = 0; i != ni_; i++) {
		for (int j = 0; j != nj_; j++) {
			pind = 9 * (lp.sizej_ * (i + 1) + j + 1);
			mind = 3 * (nj_ * i + j);

			u = data_[mind + 1];
			v = data_[mind + 2];

			func.CalculateMoment(&data_[mind], &lp.data_[pind]);

			// TODO: more options for difference definition
			diff += sqrt((u - data_[mind + 1]) * (u - data_[mind + 1]) \
				+ (v - data_[mind + 2]) * (v - data_[mind + 2]));
		}
	}

#ifdef _DEBUG

	cout << "Moment update difference: " << diff << endl;

#endif // _DEBUG
}

void LatticeMoment::OutputAscii()
{
	ofstream fout;
	string fname = "output.dat";
	fout.open(fname);

	// header
	fout << "TITLE     = \" moment \"" << endl;
	fout << "FILETYPE  = FULL" << endl;
	fout << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << ni_ << endl;
	fout << "        J = " << nj_ << endl;

	// flow variables (moments)
	int mind;
	for (int j = 0; j != nj_; j++) {
		for (int i = 0; i != ni_; i++) {
			mind = 3 * (i * nj_ + j);

			fout << left << setw(8) << i \
				<< ' ' << left << setw(8) << j \
				<< ' ' << scientific << setprecision(12) << data_[mind] \
				<< ' ' << scientific << setprecision(12) << data_[mind + 1] \
				<< ' ' << scientific << setprecision(12) << data_[mind + 2] \
				<< endl;
		}
	}

	fout.close();
}

void LatticeMoment::OutputAscii(string fname)
{
	ofstream fout;
	fout.open(fname);

	// header
	fout << "TITLE     = \" moment \"" << endl;
	fout << "FILETYPE  = FULL" << endl;
	fout << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << ni_ << endl;
	fout << "        J = " << nj_ << endl;

	// flow variables (moments)
	int mind;
	for (int j = 0; j != nj_; j++) {
		for (int i = 0; i != ni_; i++) {
			mind = 3 * (i * nj_ + j);

			fout << left << setw(8) << i \
				<< ' ' << left << setw(8) << j \
				<< ' ' << scientific << setprecision(12) << data_[mind] \
				<< ' ' << scientific << setprecision(12) << data_[mind + 1] \
				<< ' ' << scientific << setprecision(12) << data_[mind + 2] \
				<< endl;
		}
	}

	fout.close();
}
