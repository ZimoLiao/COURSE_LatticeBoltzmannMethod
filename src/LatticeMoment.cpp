#include "LatticeMoment.h"

inline int LatticeMoment::Index(int i, int j)
{
	return 3 * (nj*i + j);
}

double LatticeMoment::FuncD(double r)
{
	r = abs(r);
	if (r < 1) {
		return (3. - 2.*r + sqrt(1. + 4.*r - 4.*r*r)) / 8.;
	}
	else if (r < 2) {
		return (5. - 2.*r - sqrt(-7. + 12.*r - 4.*r*r)) / 8.;
	}
	return 0.0;
}

void LatticeMoment::Calculate(double * m, const double * f)
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

LatticeMoment::LatticeMoment()
{
	x0 = 0;
	y0 = 0;
	ni = 1;
	nj = 1;
	sizeij = 1;
	size = 3 * sizeij;

	data = new double[size];
	data[0] = 1.0;
	data[1] = 0.0;
	data[2] = 0.0;
}

LatticeMoment::~LatticeMoment()
{
	delete[] data;
}

double & LatticeMoment::operator()(int i, int j)
{
	return data[Index(i, j)];
}

void LatticeMoment::Init(int x0, int y0, int ni, int nj)
{
	this->x0 = x0;
	this->y0 = y0;
	this->ni = ni;
	this->nj = nj;
	this->sizeij = ni * nj;
	this->size = 3 * sizeij;

	data = new double[size];
	for (int i = 0; i != sizeij; i++) {
		data[3 * i + 0] = 1.0;
		data[3 * i + 1] = 0.0;
		data[3 * i + 2] = 0.0;
	}
}

void LatticeMoment::InitEntity(LatticeEntity newle)
{
	le.push_back(newle);
	nle++;
}

void LatticeMoment::Update(LatticePopulation & lp)
{
	diff = 0.0;

	double u, v;
	int mind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			mind = Index(i, j);
			u = data[mind + 1];
			v = data[mind + 2];

			Calculate(&data[mind], &lp(i, j));

			diff += sqrt((u - data[mind + 1]) * (u - data[mind + 1]) \
				+ (v - data[mind + 2]) * (v - data[mind + 2]));
		}
	}
}

void LatticeMoment::Force()
{
	for (int ie = 0; ie != nle; ie++) {

		int nm = le[ie].get_nm();
		int nm_eff = 0;

		le[ie].Reset();
		int xm, ym;
		int istart, iend, jstart, jend;

		for (int im = 0; im != nm; im++) {
			if (le[ie].IsExist(im)) {
				nm_eff++;

				xm = floor(le[ie].get_xm(im));
				ym = floor(le[ie].get_ym(im));

				istart = max(xm - 1, 0);
				iend = min(xm + 2, ni - 1);
				jstart = max(ym - 1, 0);
				jend = min(ym + 2, nj - 1);

				for (int i = istart; i <= iend; i++) {
					for (int j = jstart; j <= jend; j++) {
						le[ie].CalculateU(im, i, j, data[Index(i, j) + 1], data[Index(i, j) + 2]);
					}
				}

			}
		}

		for (int im = 0; im != nm; im++) {
			if (le[ie].IsExist(im)) {

				xm = floor(le[ie].get_xm(im));
				ym = floor(le[ie].get_ym(im));

				istart = max(xm - 1, 0);
				iend = min(xm + 2, ni - 1);
				jstart = max(ym - 1, 0);
				jend = min(ym + 2, nj - 1);

				double Fbx, Fby;
				Fbx = le[ie].CalculateFx(im, data[Index(xm, ym)]); // TODO: 这里密度取值怎么处理？
				Fby = le[ie].CalculateFy(im, data[Index(xm, ym)]);

				for (int i = istart; i <= iend; i++) {
					for (int j = jstart; j <= jend; j++) {

						data[Index(i, j) + 1] += Fbx * FuncD(double(i) - xm) \
							*FuncD(double(j) - ym) / (2.0*data[Index(i, j)]);
						data[Index(i, j) + 2] += Fby * FuncD(double(i) - xm) \
							*FuncD(double(j) - ym) / (2.0*data[Index(i, j)]);

					}
				}

			}
		}
	}
}

void LatticeMoment::WriteAscii(ofstream & fout, int step)
{
	// header
	fout << "TITLE     = Flow" << endl;
	fout << "FILETYPE  = FULL" << endl;
	fout << "VARIABLES = x y rho u v" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << ni << endl;
	fout << "        J = " << nj << endl;
	fout << "SOLUTIONTIME = " << step << endl;

	// flow variables (moments)
	int mind;
	for (int j = 0; j != nj; j++) {
		for (int i = 0; i != ni; i++) {
			mind = Index(i, j);

			fout << x0 + i << ' ' << y0 + j << ' ' \
				<< data[mind] << ' ' << data[mind + 1] << ' ' << data[mind + 2] << endl;
		}
	}
}
