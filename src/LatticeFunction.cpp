#include "LatticeFunction.h"

#include<math.h>
#include<iostream>

void LatticeFunction::CalculateEquilibrium(double* feq, double rho, double u, double v)
{
	double udotc, expr = 1.0 - 1.5 * (u * u + v * v);
	feq[0] = rho * expr * w[0];
	for (int a = 1; a != 9; a++) {
		udotc = u * cx[a] + v * cy[a];
		feq[a] = rho * (expr + udotc * (3.0 + 4.5 * udotc)) * w[a];
	}
}

void LatticeFunction::CalculateEquilibrium(double* feq, double* moment)
{

#ifdef _DEBUG
	if (fabs(moment[0] - 1.0) > 1.0) {
		std::cout << "out_of_range\n";
		return;
	}
#endif // _DEBUG

	double udotc, expr = 1.0 - 1.5 * (moment[1] * moment[1] + moment[2] * moment[2]);
	feq[0] = moment[0] * expr * w[0];
	for (int a = 1; a != 9; a++) {
		udotc = moment[1] * cx[a] + moment[2] * cy[a];
		feq[a] = moment[0] * (expr + udotc * (3.0 + 4.5 * udotc)) * w[a];
	}
}
