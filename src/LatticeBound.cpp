#include "LatticeBound.h"

<<<<<<< HEAD
LatticeBound::LatticeBound(int index_global, int type, int direct, int is, int ie, int js, int je, double rho, double u, double v)
{
	this->index_global = index_global;
	this->type = type;
	this->direct = direct;
	this->is = is;
	this->ie = ie;
	this->js = js;
	this->je = je;

	this->rho = rho;
	this->u = u;
	this->v = v;
}

void LatticeBound::UpdateLocal(int rank, int ni)
{
	int isr = ni * rank, ier = (ni + 1)*rank - 1;
	is = std::max(is, isr);
	ie = std::max(ie, ier);
}

bool LatticeBound::IsExist(int rank, int ni)
{
	int isr = ni * rank, ier = (ni + 1)*rank - 1;

	return !(is > ier || ie < isr);
=======
LatticeBound::LatticeBound(int rank, int ni, int nj, \
	int type, int direct, \
	int istart, int iend, int jstart, int jend, \
	double rho, double u, double v)
{
	imin = rank * ni;
	imax = (rank + 1)*ni - 1;
	jmin = 0;
	jmax = nj;

	this->type = type;
	this->direct = direct;
	this->istart = istart;
	this->iend = iend;
	this->jstart = jstart;
	this->jend = jend;
	this->rho = rho;
	this->u = u;
	this->v = v;

	exist = !(jend<jmin || iend<imin || jstart>jmax || istart>imax);

	this->istart = max(istart, imin) % ni;
	this->iend = min(iend, imax) % ni;
	this->jstart = max(jstart, jmin) % nj;
	this->jend = min(jend, jmax) % nj;
}

bool LatticeBound::IsExist()
{
	return exist;
>>>>>>> new
}

void LatticeBound::CalculateNebb0(double * f)
{
	switch (direct)
	{
	case 1:
		f[1] = f[3];
		f[5] = f[7] - 0.5*(f[2] - f[4]);
		f[8] = f[6] + 0.5*(f[2] - f[4]);
		break;
	case 3:
		f[3] = f[1];
		f[6] = f[8] - 0.5*(f[2] - f[4]);
		f[7] = f[5] + 0.5*(f[2] - f[4]);
		break;
	case 2:
		f[2] = f[4];
		f[5] = f[7] - 0.5*(f[1] - f[3]);
		f[6] = f[8] + 0.5*(f[1] - f[3]);
		break;
	case 4:
		f[4] = f[2];
		f[8] = f[6] - 0.5*(f[1] - f[3]);
		f[7] = f[5] + 0.5*(f[1] - f[3]);
		break;
	case 5:
		f[1] = f[3];
		f[2] = f[4];
		f[5] = f[7];
		f[6] = 0.5*(rho - f[0] - 2.0*(f[3] + f[4] + f[7]));
		f[8] = f[6];
		break;
	case 6:
		f[3] = f[1];
		f[2] = f[4];
		f[6] = f[8];
		f[5] = 0.5*(rho - f[0] - 2.0*(f[1] + f[4] + f[8]));
		f[7] = f[5];
		break;
	case 7:
		f[3] = f[1];
		f[4] = f[2];
		f[7] = f[5];
		f[6] = 0.5*(rho - f[0] - 2.0*(f[1] + f[2] + f[5]));
		f[8] = f[6];
		break;
	case 8:
		f[1] = f[3];
		f[4] = f[2];
		f[8] = f[6];
		f[5] = 0.5*(rho - f[0] - 2.0*(f[2] + f[3] + f[6]));
		f[7] = f[5];
		break;
	}
}

<<<<<<< HEAD
=======
void LatticeBound::CalculateNebb0(double * f, double rho_ref)
{
	switch (direct)
	{
	case 1:
		f[1] = f[3];
		f[5] = f[7] - 0.5*(f[2] - f[4]);
		f[8] = f[6] + 0.5*(f[2] - f[4]);
		break;
	case 3:
		f[3] = f[1];
		f[6] = f[8] - 0.5*(f[2] - f[4]);
		f[7] = f[5] + 0.5*(f[2] - f[4]);
		break;
	case 2:
		f[2] = f[4];
		f[5] = f[7] - 0.5*(f[1] - f[3]);
		f[6] = f[8] + 0.5*(f[1] - f[3]);
		break;
	case 4:
		f[4] = f[2];
		f[8] = f[6] - 0.5*(f[1] - f[3]);
		f[7] = f[5] + 0.5*(f[1] - f[3]);
		break;
	case 5:
		f[1] = f[3];
		f[2] = f[4];
		f[5] = f[7];
		f[6] = 0.5*(rho_ref - f[0] - 2.0*(f[3] + f[4] + f[7]));
		f[8] = f[6];
		break;
	case 6:
		f[3] = f[1];
		f[2] = f[4];
		f[6] = f[8];
		f[5] = 0.5*(rho_ref - f[0] - 2.0*(f[1] + f[4] + f[8]));
		f[7] = f[5];
		break;
	case 7:
		f[3] = f[1];
		f[4] = f[2];
		f[7] = f[5];
		f[6] = 0.5*(rho_ref - f[0] - 2.0*(f[1] + f[2] + f[5]));
		f[8] = f[6];
		break;
	case 8:
		f[1] = f[3];
		f[4] = f[2];
		f[8] = f[6];
		f[5] = 0.5*(rho_ref - f[0] - 2.0*(f[2] + f[3] + f[6]));
		f[7] = f[5];
		break;
	}
}

>>>>>>> new
void LatticeBound::CalculateNebbV(double * f)
{
	switch (direct)
	{
	case 1:
		rho = (f[0] + f[2] + f[4] + 2.0*(f[3] + f[6] + f[7])) / (1.0 - u);

		f[1] = f[3] + 2.0*rho*u / 3.0;
		f[5] = f[7] - 0.5*(f[2] - f[4]) + 0.5*rho*v + rho * u / 6.0;
		f[8] = f[6] + 0.5*(f[2] - f[4]) - 0.5*rho*v + rho * u / 6.0;
		break;
	case 3:
		rho = (f[0] + f[2] + f[4] + 2.0*(f[1] + f[5] + f[8])) / (1.0 - u);

		f[3] = f[1] - 2.0*rho*u / 3.0;
		f[6] = f[8] - 0.5*(f[2] - f[4]) + 0.5*rho*v - rho * u / 6.0;
		f[7] = f[5] + 0.5*(f[2] - f[4]) - 0.5*rho*v - rho * u / 6.0;
		break;
	case 2:
		rho = (f[0] + f[1] + f[3] + 2.0*(f[4] + f[7] + f[8])) / (1.0 - v);

		f[2] = f[4] + 2.0*rho*v / 3.0;
		f[5] = f[7] - 0.5*(f[1] - f[3]) + 0.5*rho*u + rho * v / 6.0;
		f[6] = f[8] + 0.5*(f[1] - f[3]) - 0.5*rho*u + rho * v / 6.0;
		break;
	case 4:
		rho = (f[0] + f[1] + f[3] + 2.0*(f[2] + f[5] + f[6])) / (1.0 - v);

		f[4] = f[2] - 2.0*rho*v / 3.0;
		f[8] = f[6] - 0.5*(f[1] - f[3]) + 0.5*rho*u - rho * v / 6.0;
		f[7] = f[5] + 0.5*(f[1] - f[3]) - 0.5*rho*u - rho * v / 6.0;
		break;
	}
}

void LatticeBound::CalculateNebbP(double * f)
{
	double jx, jy;
	switch (direct)
	{
	case 1:
		jx = rho - (f[0] + f[2] + f[4] + 2.0*(f[3] + f[6] + f[7]));

		f[1] = f[3] + 2.0*jx / 3.0;
		f[5] = f[7] - 0.5*(f[2] - f[4]) + jx / 6.0;
		f[8] = f[6] + 0.5*(f[2] - f[4]) + jx / 6.0;
		break;
	case 3:
		jx = (f[0] + f[2] + f[4] + 2.0*(f[1] + f[5] + f[8])) - rho;

		f[3] = f[1] - 2.0*jx / 3.0;
		f[6] = f[8] - 0.5*(f[2] - f[4]) - jx / 6.0;
		f[7] = f[5] + 0.5*(f[2] - f[4]) - jx / 6.0;
		break;
	case 2:
		jy = rho - (f[0] + f[1] + f[3] + 2.0*(f[4] + f[7] + f[8]));

		f[2] = f[4] + 2.0*jy / 3.0;
		f[5] = f[7] - 0.5*(f[1] - f[3]) + jy / 6.0;
		f[6] = f[8] + 0.5*(f[1] - f[3]) + jy / 6.0;
		break;
	case 4:
		jy = (f[0] + f[1] + f[3] + 2.0*(f[2] + f[5] + f[6])) - rho;

		f[4] = f[2] - 2.0*jy / 3.0;
		f[8] = f[6] - 0.5*(f[1] - f[3]) - jy / 6.0;
		f[7] = f[5] + 0.5*(f[1] - f[3]) - jy / 6.0;
		break;
	}
}
<<<<<<< HEAD
=======

int LatticeBound::GetType()
{
	return type;
}

int LatticeBound::GetIs()
{
	return istart;
}

int LatticeBound::GetIe()
{
	return iend;
}

int LatticeBound::GetJs()
{
	return jstart;
}

int LatticeBound::GetJe()
{
	return jend;
}
>>>>>>> new
