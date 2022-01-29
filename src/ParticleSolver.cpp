#include "ParticleSolver.h"

void ParticleSolver::Ignore(ifstream& fin, int l)
{
	char buffer[101];
	for (int i = 0; i != l; i++) {
		fin.getline(buffer, 100);
	}
}

ParticleSolver::ParticleSolver()
{
	Force = new double[1];
	buffer_host = new double[1];
}

ParticleSolver::~ParticleSolver()
{
	delete[] Force;
	fout.close();// TODO 不应该这样偷懒。。。需修改
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

		ImmersedParticle newpart(rank * ni, 0, ni, nj, r, x0, y0, phi0, ux0, uy0, uphi0);
		part.push_back(newpart);
	}

	fin.close();

	/* force calculation */
	MPI_Comm_size(MPI_COMM_WORLD, &size); // TODO
	if (rank == host) {
		buffer_host = new double[size * npart * 3];// TODO 没有考虑各颗粒尺寸不同的情况
	}
	Force = new double[npart * 3];

	/* write to file */
	fout.open("out/force.dat");
	fout << "variables = step Fx Fy M\n";// 只适用于一个颗粒
}

void ParticleSolver::Calculate(LatticeMoment& lm)
{
	for (int ip = 0; ip != npart; ip++) {
		for (int iforce = 0; iforce != nforce; iforce++) {
			if (part[ip].IsExist()) {
				part[ip].ForceMoment(lm);
			}
			lm.UpdateGhost();
		}
		part[ip].CalculateTotalForce();
	}
}

void ParticleSolver::UpdateForce()
{
	step++;
	for (int i = 0; i != npart; i++) {
		Force[3 * i] = part[i].Fx;
		Force[3 * i + 1] = part[i].Fy;
		Force[3 * i + 2] = part[i].M;
	}
	MPI_Gather(Force, npart * 3, MPI_DOUBLE, buffer_host, npart * 3, MPI_DOUBLE, host, MPI_COMM_WORLD);
	if (rank == host) {
		fout << step << ' ';
		for (int i = 0; i != npart * 3; i++) {
			Force[i] = 0.0;
			for (int j = 0; j != size; j++) {
				Force[i] += buffer_host[i + j * 3 * npart];
			}
			fout << std::scientific << std::setprecision(12) << Force[i] << '\t';
		}
		fout << '\n';
	}
}
