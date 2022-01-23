#include "LatticeEntity.h"

void LatticeEntity::UpdateMarker()
{
	for (int i = 0; i != nm; i++) {
		xm[i] = x + r * cos(dtheta*i);
		ym[i] = x + r * sin(dtheta*i);
	}
}

LatticeEntity::LatticeEntity(double r, double x, double y)
{
	this->r = r;
	nm = floor(2 * pi*this->r);
	dtheta = theta / double(nm);

	xm = new double[nm];
	ym = new double[nm];

	theta = 0.0;
	this->x = x;
	this->y = y;
}

LatticeEntity::~LatticeEntity()
{
	delete[] xm;
	delete[] ym;
}
