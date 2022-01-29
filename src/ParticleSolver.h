#ifndef PARTICLESOLVER_H_
#define PARTICLESOLVER_H_

#include<fstream>
#include<string>

#include"ImmersedParticle.h"
#include"LatticeMoment.h"

using namespace std;

class ParticleSolver
{
	const int host = 0;

	int rank, size;
	int ni, nj;

	int step = 0;
	double* Force;
	double* buffer_host;


	/* file stream (write force) */
	ofstream fout;


	/* internal functions */
	void Ignore(ifstream& fin, int l);

public:
	int nforce = 0;
	int npart = 0;
	vector<ImmersedParticle> part;

	/* constructor */
	ParticleSolver();
	~ParticleSolver();


	void Init(int rank, int ni, int nj);

	void Calculate(LatticeMoment& lm);

	void UpdateForce();


};


#endif // !PARTICLESOLVER_H_