#include "LatticeFunction.h"

#include<math.h>
#include<iostream>

void LatticeFunction::CalculateEquilibrium(double* feq, const double rho, const double u, const double v)
{
	double udotc, expr = 1.0 - 1.5 * (u * u + v * v);
	feq[0] = rho * expr * w[0];
	for (int a = 1; a != 9; a++) {
		udotc = u * cx[a] + v * cy[a];
		feq[a] = rho * (expr + udotc * (3.0 + 4.5 * udotc)) * w[a];
	}
}

void LatticeFunction::CalculateEquilibrium(double* feq, const double* m)
{

#ifdef _DEBUG
	if (fabs(m[0] - 1.0) > 1.0) {
		std::cout << "out_of_range\n";
		return;
	}
#endif // _DEBUG

	double udotc, expr = 1.0 - 1.5 * (m[1] * m[1] + m[2] * m[2]);
	feq[0] = m[0] * expr * w[0];
	for (int a = 1; a != 9; a++) {
		udotc = m[1] * cx[a] + m[2] * cy[a];
		feq[a] = m[0] * (expr + udotc * (3.0 + 4.5 * udotc)) * w[a];
	}
}

void LatticeFunction::CalculateMoment(double* m, const double* f)
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
