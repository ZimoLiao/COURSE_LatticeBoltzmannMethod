#include "LatticePopulation.h"

LatticePopulation::~LatticePopulation()
{
	delete[] data;
}

void LatticePopulation::InitBoundary(LatticeBound newlb)
{
	lb.push_back(newlb);
	nlb++;
}
