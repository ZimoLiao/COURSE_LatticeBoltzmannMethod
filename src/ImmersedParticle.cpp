#include "ImmersedParticle.h"

inline bool ImmersedParticle::IsExist(double mxi)
{
	return !(mxi <= xmin || mxi >= xmax);
}

double ImmersedParticle::D(double d)
{
	d = abs(d);
	if (d < 1.0) {
		return (3. - 2. * d + sqrt(1. + 4. * d - 4. * d * d)) / 8.;
	}
	else if (d < 2.0) {
		return (5. - 2. * d - sqrt(-7. + 12. * d - 4. * d * d)) / 8.;
	}
	return 0.0;
}

double ImmersedParticle::Delta(int im, double i, double j)
{
	return D(mx[im] - i) * D(my[im] - j);
}

ImmersedParticle::ImmersedParticle(\
	double x0, double y0, double ni, double nj, double r, \
	double gx, double gy, double gphi, \
	double ux, double uy, double uphi)
{
	this->r = r;
	this->x = gx - x0;
	this->y = gy - y0;
	this->phi = gphi;
	this->ux = ux;
	this->uy = uy;
	this->uphi = uphi;

	xmin = -2.0;
	ymin = -2.0;
	xmax = ni + 1.0;
	ymax = nj + 1.0;

	/* markers initialization */
	nm = ceil(2.0 * pi * r);
	dphi = 2.0 * pi / double(nm);

	// allocation
	mx = new double[nm];
	my = new double[nm];
	mux = new double[nm];
	muy = new double[nm];
	mufx = new double[nm];
	mufy = new double[nm];
	mFx = new double[nm];
	mFy = new double[nm];
	mexist = new bool[nm];

	// value calculation
	double phii;
	exist = false;
	for (int i = 0; i != nm; i++) {
		phii = phi + dphi * i;

		mx[i] = x + r * cos(phii);
		my[i] = y + r * sin(phii);
		mux[i] = ux - uphi * r * sin(phii);
		muy[i] = uy + uphi * r * cos(phii);
		mufx[i] = 0.0;
		mufy[i] = 0.0;
		mFx[i] = 0.0;
		mFy[i] = 0.0;

		mexist[i] = IsExist(mx[i]);
		if (mexist[i]) { exist = true; }
	}
}

ImmersedParticle::ImmersedParticle(const ImmersedParticle& newpart)
{
	this->r = newpart.r;
	this->x = newpart.x;
	this->y = newpart.y;
	this->phi = newpart.phi;
	this->ux = newpart.ux;
	this->uy = newpart.uy;
	this->uphi = newpart.uphi;

	xmin = newpart.xmin;
	ymin = newpart.ymin;
	xmax = newpart.xmax;
	ymax = newpart.ymax;

	/* markers initialization */
	nm = newpart.nm;
	dphi = newpart.dphi;

	// allocation
	mx = new double[nm];
	my = new double[nm];
	mux = new double[nm];
	muy = new double[nm];
	mufx = new double[nm];
	mufy = new double[nm];
	mFx = new double[nm];
	mFy = new double[nm];
	mexist = new bool[nm];

	// value calculation
	exist = newpart.exist;
	for (int i = 0; i != nm; i++) {
		mx[i] = newpart.mx[i];
		my[i] = newpart.my[i];
		mux[i] = newpart.mux[i];
		muy[i] = newpart.muy[i];
		mufx[i] = newpart.mufx[i];
		mufy[i] = newpart.mufy[i];
		mFx[i] = newpart.mFx[i];
		mFy[i] = newpart.mFy[i];
		mexist[i] = newpart.mexist[i];
	}
}

ImmersedParticle::~ImmersedParticle()
{
	delete[] mx;
	delete[] my;
	delete[] mux;
	delete[] muy;
	delete[] mufx;
	delete[] mufy;
	delete[] mFx;
	delete[] mFy;
	delete[] mexist;
}

bool ImmersedParticle::IsExist()
{
	return exist;
}

void ImmersedParticle::ForceMoment(LatticeMoment& lm)
{
	/* update unforced velocity */
	int istart, iend, jstart, jend;

	for (int im = 0; im != nm; im++) {
		mufx[im] = 0.;
		mufy[im] = 0.;

		if (mexist[im]) {

			istart = floor(mx[im]) - 1;
			jstart = floor(my[im]) - 1;
			iend = istart + 3;
			jend = jstart + 3;

			double delta = 0.0;
			for (int i = istart; i != iend; i++) {
				for (int j = jstart; j != jend; j++) {
					delta = Delta(im, i, j);
					mufx[im] += delta * lm(i, j, 1);
					mufy[im] += delta * lm(i, j, 2);
				}
			}

		}
	}

	/* update force and moment */

	for (int im = 0; im != nm; im++) {
		mFx[im] = 0.;
		mFy[im] = 0.;

		if (mexist[im]) {

			// density interpolation
			double rho = Interpolation(im, lm);

			mFx[im] = 2.0 * rho * (mux[im] - mufx[im]);
			mFy[im] = 2.0 * rho * (muy[im] - mufy[im]);

			istart = floor(mx[im]) - 1;
			jstart = floor(my[im]) - 1;
			iend = istart + 3;
			jend = jstart + 3;

			double delta = 0.0;
			for (int i = istart; i != iend; i++) {
				for (int j = jstart; j != jend; j++) {
					delta = Delta(im, i, j);
					rho = lm(i, j, 0);

					lm(i, j, 1) += delta * mFx[im] / 2.0 / rho;
					lm(i, j, 2) += delta * mFy[im] / 2.0 / rho;
				}
			}

		}
	}
}

void ImmersedParticle::CalculateTotalForce()
{
	Fx = 0.0;
	Fy = 0.0;
	M = 0.0;

	double ds = r * dphi;
	if (IsExist()) {
		for (int im = 0; im != nm; im++) {
			if (mexist[im]) {
				if (mx[im] < xmax - 2.0 && mx[im]>0.0) {
					Fx -= mFx[im] * ds;
					Fy -= mFy[im] * ds;
					M += -mFy[im] * (mx[im] - x) + mFx[im] * (my[im] - y) * ds;
				}
			}
		}
	}
}

double ImmersedParticle::Interpolation(int im, LatticeMoment& lm)
{
	double dx, dy, dxc, dyc;
	int i0, i1, j0, j1;

	i0 = floor(mx[im]);
	j0 = floor(my[im]);
	i1 = ceil(mx[im]);
	j1 = ceil(my[im]);
	dx = mx[im] - i0;
	dxc = 1.0 - dx;
	dy = my[im] - j0;
	dyc = 1.0 - dy;

	return lm(i0, j0, 0) * dxc * dyc + \
		lm(i1, j0, 0) * dx * dyc + \
		lm(i0, j1, 0) * dxc * dy + \
		lm(i1, j1, 0) * dx * dy;
}

int ImmersedParticle::GetNm()
{
	return nm;
}
