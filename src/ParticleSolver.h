#ifndef PARTICLESOLVER_H_
#define PARTICLESOLVER_H_

#include<fstream>
#include<string>

#include"ImmersedParticle.h"

using namespace std;

class ParticleSolver
{
	int rank;
	int ni, nj;

	/* internal functions */
	void Ignore(ifstream& fin, int l);

public:
	int nforce = 0;
	int npart = 0;
	vector<ImmersedParticle> part;

	/* constructor */
	ParticleSolver();


	void Init(int rank, int ni, int nj);

};


#endif // !PARTICLESOLVER_H_