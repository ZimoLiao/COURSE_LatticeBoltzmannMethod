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

public: // TODO: 斟酌参数的可见性

	// collision operator
	// SRT
	double tau = 0.54, somega = 1.0 / tau, somegac = 1.0 - somega;

	// MRT (with Gram-Schmidt orthogonal moment)
	// TODO: 检查是否有误
	double momega1 = 0.03, momega2 = 1.0, momega3 = 1.0, momega4 = 1.97;
	double momega[9] = { 0., momega1, momega2, 0., momega3, 0., momega3, momega4, momega4 };
	double momegac[9] = { 1., 1. - momega1, 1. - momega2, 1., 1. - momega3, \
		1., 1. - momega3, 1. - momega4, 1. - momega4 };
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

	/* functions */
	void CalculateFeq(double* feq, const double rho, const double u, const double v);
	void CalculateFeq(double* feq, const double* m);
	void CalculateFstar(double* f, const double* m); // MRT
	void CalculateMoment(double* m, const double* f);

	void Init(double tau, double momega1, double momega2, double momega3, double momega4);
};


#endif // !LATTICEFUNCTION_H_