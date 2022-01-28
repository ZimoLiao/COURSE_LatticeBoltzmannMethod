#include "ParticleSolver.h"

void ParticleSolver::Ignore(ifstream & fin, int l)
{
	char buffer[101];
	for (int i = 0; i != l; i++) {
		fin.getline(buffer, 100);
	}
}

ParticleSolver::ParticleSolver()
{
}

void ParticleSolver::Init(int rank, int ni, int nj)
{
	this->rank = rank;
	this->ni = ni;
	this->nj = nj;

	ifstream fin;
	fin.open("in/part.init");

	Ignore(fin, 4);

	fin >> nforce;
	Ignore(fin, 4);

	fin >> npart;
	Ignore(fin, 2);

	double r, x0, y0, phi0, ux0, uy0, uphi0;
	for (int i = 0; i != npart; i++) {
		fin >> r >> x0 >> y0 >> phi0 >> ux0 >> uy0 >> uphi0;
		Ignore(fin, 1);

		ImmersedParticle newpart(rank*ni, 0, ni, nj, r, x0, y0, phi0, ux0, uy0, uphi0);
		part.push_back(newpart);
	}

	fin.close();
}
