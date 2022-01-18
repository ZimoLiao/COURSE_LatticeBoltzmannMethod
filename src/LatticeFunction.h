#ifndef LATTICEFUNCTION_H_
#define LATTICEFUNCTION_H_

class LatticeFunction
{
	/* parameters */

	// velocity set
	const double cx[9] = { 0.,1.,0.,-1.,0.,1.,-1.,-1.,1. };
	const double cy[9] = { 0.,0.,1.,0.,-1.,1.,1.,-1.,-1. };
	const double w[9] = {
		4. / 9.,
		1. / 9., 1. / 9., 1. / 9., 1. / 9.,
		1. / 36., 1. / 36., 1. / 36., 1. / 36.
	};

	// collision operator
	// SRT
	double tau = 0.8, omega = 1.0 / tau, omegac = 1.0 - omega;

	// MRT (with Gram-Schmidt orthogonal moment)
	// TODO: ºÏ≤È «∑Ò”–ŒÛ
	double omega1 = 1.0, omega2 = 1.0, omega3 = 1.0, omega4 = 1.0;
	const double M[9][9] = {
		{	1,	1,	1,	1,	1,	1,	1,	1,	1	},
		{	-4,	-1,	-1,	-1,	-1,	2,	2,	2,	2	},
		{	4,	-2,	-2,	-2,	-2,	1,	1,	1,	1	},
		{	0,	1,	0,	-1,	0,	1,	-1,	-1,	1	},
		{	0,	-2,	0,	2,	0,	1,	-1,	-1,	1	},
		{	0,	0,	1,	0,	-1,	1,	1,	-1,	-1	},
		{	0,	0,	-2,	0,	2,	1,	1,	-1,	-1	},
		{	0,	1,	-1,	1,	-1,	0,	0,	0,	0	},
		{	0,	0,	0,	0,	0,	1,	-1,	1,	-1	},
	};
	const double Minv[9][9] = {
		{	1. / 9,	-1. / 9,	1. / 9,		0,			0,			0,			0,			0,			0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	1. / 6,		-1. / 6,	0,			0,			1. / 4,		0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	0,			0,			1. / 6,		-1. / 6,	-1. / 4,	0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	-1. / 6,	1. / 6,		0,			0,			1. / 4,		0		},
		{	1. / 9,	-1. / 36,	-1. / 18,	0,			0,			-1. / 6,	1. / 6,		-1. / 4,	0		},
		{	1. / 9,	1. / 18,	1. / 36,	1. / 6,		1. / 12,	1. / 6,		1. / 12,	0,			1. / 4	},
		{	1. / 9,	1. / 18,	1. / 36,	-1. / 6,	-1. / 12,	1. / 6,		1. / 12,	0,			-1. / 4	},
		{	1. / 9,	1. / 18,	1. / 36,	-1. / 6,	-1. / 12,	-1. / 6,	-1. / 12,	0,			1. / 4	},
		{	1. / 9,	1. / 18,	1. / 36,	1. / 6,		1. / 12,	-1. / 6,	-1. / 12,	0,			-1. / 4	},
	};

public:


	/* functions */
	void CalculateEquilibrium(double* feq, const double rho, const double u, const double v);
	void CalculateEquilibrium(double* feq, const double* m);
	void CalculateMoment(double* m, const double* f);
};


#endif // !LATTICEFUNCTION_H_