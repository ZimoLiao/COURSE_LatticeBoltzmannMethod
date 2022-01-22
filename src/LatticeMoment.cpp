#include "LatticeMoment.h"

void LatticeMoment::CalculateMoment(double * m, const double * f)
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

inline int LatticeMoment::IndexM(int i, int j)
{
	return 3 * (nj * i + j);
}

LatticeMoment::LatticeMoment()
{
	i0 = 0;
	j0 = 0;
	ni = 1;
	nj = 1;
	sizeij = 1;
	size = 3 * sizeij;

	data = new double[size];
	data[0] = 1.0;
	data[1] = 0.0;
	data[2] = 0.0;
}

LatticeMoment::LatticeMoment(const LatticeMoment & lm)
{
	i0 = lm.i0;
	j0 = lm.j0;
	ni = lm.ni;
	nj = lm.nj;
	sizeij = lm.sizeij;
	size = lm.size;

	data = new double[size];
	for (int i = 0; i != size; i++) { data[i] = lm.data[i]; }
}

LatticeMoment::~LatticeMoment()
{
	delete[] data;
}

void LatticeMoment::Init(int i0, int j0, int ni, int nj)
{
	this->i0 = i0;
	this->j0 = j0;
	this->ni = ni;
	this->nj = nj;
	this->sizeij = ni * nj;
	this->size = 3 * sizeij;

	data = new double[size];
	for (int i = 0; i != sizeij; i++) {
		data[i + 0] = 1.0;
		data[i + 1] = 0.0;
		data[i + 2] = 0.0;
	}
}

// TODO
void LatticeMoment::Update(LatticePopulation & lp)
{
	diff = 0.0;

	int find, mind;
	double u, v;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			mind = IndexM(i, j);
			find = lp.IndexF(i + 1, j + 1);

			u = data[mind + 1];
			v = data[mind + 2];

			CalculateMoment(&data[mind], &lp.data[find]);

			diff += sqrt((u - data[mind + 1]) * (u - data[mind + 1]) \
				+ (v - data[mind + 2]) * (v - data[mind + 2]));
		}
	}
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
	fout << "        I = " << ni << endl;
	fout << "        J = " << nj << endl;
	fout << "SOLUTIONTIME = " << step << endl;

	// flow variables (moments)
	int mind;
	for (int j = 0; j != nj; j++) {
		for (int i = 0; i != ni; i++) {
			mind = IndexM(i, j);

			fout << std::left << std::setw(8) << i0 + i \
				<< ' ' << std::left << std::setw(8) << j0 + j \
				<< ' ' << std::scientific << std::setprecision(12) << data[mind] \
				<< ' ' << std::scientific << std::setprecision(12) << data[mind + 1] \
				<< ' ' << std::scientific << std::setprecision(12) << data[mind + 2] \
				<< endl;
		}
	}

	fout.close();
}
