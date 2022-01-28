#ifndef LATTICEMOMENT_H_
#define LATTICEMOMENT_H_

#include<string>
#include<fstream>
#include<iomanip>

#include"LatticeBlock.h"
#include"LatticePopulation.h"
#include"ImmersedParticle.h"

using std::string;

constexpr double pi = 3.1415926535897932;

class LatticeMoment :
	public LatticeBlock<3, 3>
{



public:
	friend class LatticePopulation;
	friend class ImmersedParticle;


	/* monitors */
	double diff = 1.0;


	/* initialization */
	void InitData();

	// TODO: just for test
	void InitDataShear();
	// TODO: need delete


	/* functions */
	void Update(LatticePopulation& lp);

	void WriteAscii(int step);

};


#endif // !LATTICEMOMENT_H_