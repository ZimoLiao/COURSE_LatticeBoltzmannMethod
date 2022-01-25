#include "LatticeMoment.h"

LatticeMoment::~LatticeMoment()
{
	delete[] data;
}

void LatticeMoment::InitEntity(LatticeEntity newle)
{
	le.push_back(newle);
	nle++;
}
