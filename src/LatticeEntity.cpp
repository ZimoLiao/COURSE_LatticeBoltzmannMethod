#include "LatticeEntity.h"

inline bool LatticeEntity::IsExist(int im)
{
	return xm[im] > -2.0 && xm[im] < (xe - xs + 1.0);
}

void LatticeEntity::CalculateU(int im, double i, double j, double u, double v)
{
	double delta = FuncD(xm[im] - i)*FuncD(ym[im] - j);
	uxf[im] += delta * u;
	uyf[im] += delta * v;
}

double LatticeEntity::CalculateFx(int im, double rho)
{
	return 2.0*rho*(uxm[im] - uxf[im]);
}

double LatticeEntity::CalculateFy(int im, double rho)
{
	return 2.0*rho*(uym[im] - uyf[im]);
}

void LatticeEntity::Reset()
{
	for (int i = 0; i != nm; i++) {
		uxf[i] = 0.0;
		uyf[i] = 0.0;
	}
}

double LatticeEntity::FuncD(double r)
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

LatticeEntity::LatticeEntity(int index_global, double r, double xg, double yg, double phi, double xs, double xe, double ys, double ye, double ux0, double uy0, double uphi0)
{
	this->index_global = index_global;

	this->r = r;
	this->xg = xg;
	this->yg = yg;
	this->phi = phi;

	this->xs = xs;
	this->xe = xe;
	this->ys = ys;
	this->ye = ye;

	x = xg - xs;
	y = yg;
	ux = ux0;
	uy = uy0;
	uphi = uphi0;

	nm = ceil(2.0*pi*r);
	dphi = 2.0*pi / double(nm);

	xm = new double[nm];
	ym = new double[nm];
	uxm = new double[nm];
	uym = new double[nm];
	uxf = new double[nm];
	uyf = new double[nm];

	for (int i = 0; i != nm; i++) {
		xm[i] = x + r * cos(dphi*i + phi);
		ym[i] = y + r * sin(dphi*i + phi);
		uxm[i] = ux - uphi * r * sin(dphi*i + phi);
		uym[i] = uy + uphi * r * cos(dphi*i + phi);
		uxf[i] = 0.0;
		uyf[i] = 0.0;
	}
}

LatticeEntity::LatticeEntity(const LatticeEntity & le)
{
	this->index_global = le.index_global;

	this->r = le.r;
	this->xg = le.xg;
	this->yg = le.yg;
	this->phi = le.phi;

	this->xs = le.xs;
	this->xe = le.xe;
	this->ys = le.ys;
	this->ye = le.ye;

	x = le.x;
	y = le.y;
	ux = le.ux;
	uy = le.uy;
	uphi = le.uphi;

	nm = le.nm;
	dphi = le.dphi;

	xm = new double[nm];
	ym = new double[nm];
	uxm = new double[nm];
	uym = new double[nm];
	uxf = new double[nm];
	uyf = new double[nm];

	for (int i = 0; i != nm; i++) {
		xm[i] = le.xm[i];
		ym[i] = le.ym[i];
		uxm[i] = le.uxm[i];
		uym[i] = le.uym[i];
		uxf[i] = 0.0;
		uyf[i] = 0.0;
	}
}

LatticeEntity::~LatticeEntity()
{
	delete[] xm;
	delete[] ym;
	delete[] uxm;
	delete[] uym;
	delete[] uxf;
	delete[] uyf;
}

bool LatticeEntity::IsExist()
{
	for (int i = 0; i != nm; i++) {
		if (IsExist(i)) {
			return true;
		}
	}
	return false;
}

int LatticeEntity::get_nm()
{
	return nm;
}

double LatticeEntity::get_xm(int im)
{
	return xm[im];
}

double LatticeEntity::get_ym(int im)
{
	return ym[im];
}

double LatticeEntity::get_uxm(int im)
{
	return uxm[im];
}

double LatticeEntity::get_uym(int im)
{
	return uym[im];
}
