#include "LatticeSolver.h"

void LatticeSolver::Ignore(int l)
{
	char buffer[101];
	for (int i = 0; i != l; i++) {
		fin.getline(buffer, 100);
	}
}
LatticeSolver::LatticeSolver()
{
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == host) {
		// read input
		fin.open("in/input");

		Ignore(4);
		fin >> ni >> nj;
		Ignore(3);
		fin >> di >> dj;
		Ignore(3);
		for (int i = 1; i != 9; i++) {
			fin >> b[i];
		}
		Ignore(5);
		fin >> is_srt >> tau;
		Ignore(2);
		fin >> is_mrt >> omega_e >> omega_ep >> omega_q >> omega_nu;

		fin.close();
	}

	// broadcast parameters 
	int buffer1[14] = { ni, nj,di, dj, b[1], b[2] ,b[3] ,b[4] ,b[5] ,b[6] ,b[7] ,b[8] ,is_srt, is_mrt };
	double buffer2[5] = { tau, omega_e, omega_ep, omega_q, omega_nu };
	MPI_Bcast(buffer1, 14, MPI_INT, host, MPI_COMM_WORLD);
	MPI_Bcast(buffer2, 5, MPI_DOUBLE, host, MPI_COMM_WORLD);

	if (rank != host) {
		ni = buffer1[0];
		nj = buffer1[1];
		di = buffer1[2];
		dj = buffer1[3];
		for (int i = 1; i != 9; i++) {
			b[i] = buffer1[i + 3];
		}

		is_srt = buffer1[12];
		is_mrt = buffer1[13];

		tau = buffer2[0];
		omega_e = buffer2[1];
		omega_ep = buffer2[2];
		omega_q = buffer2[3];
		omega_nu = buffer2[4];
	}

	// populations and moments initialization
	// dj==1
	lm.Init(rank*ni / di, 0, ni / di, nj / dj);
	lp.Init(lm);

	lp.InitParameter(tau);
	lp.InitConnection(di, b);

	// read in boundaries
	fin.open("in/input");

	Ignore(22);
	int nlb, nlbrank = 0;
	fin >> nlb;
	Ignore(2);
	for (int b = 0; b != nlb; b++) {
		int type, direct, is, ie, js, je;
		double rho, u, v;

		fin >> type >> direct >> is >> ie >> js >> je >> rho >> u >> v;
		Ignore(1);

		int isrank = rank * ni / di, ierank = (rank + 1)*ni / di, jsrank, jerank;

		if (!(ie <= isrank || is >= ierank)) {
			isrank = max(is, isrank);
			ierank = min(ie, ierank);

			isrank = isrank % (ni / di) + 1;
			ierank = ierank % (ni / di);
			if (ierank == 0) { ierank += ni / di; }
			ierank++;
			jsrank = js + 1;
			jerank = je + 1;

			LatticeBound lb(type, direct, isrank, ierank, jsrank, jerank, rho, u, v);
			lp.InitBoundary(lb);

			nlbrank++;
		}
	}
	cout << "rank " << rank << " has " << nlbrank << " boundary units.\n";

	fin.close();
}

void LatticeSolver::Calculate()
{
	if (is_srt) {
		lp.CollideSrt(lm);
	}
	lp.UpdateGhost();
	lp.Stream();
	lp.Boundary();

	lm.Update(lp);
}

void LatticeSolver::OutputAscii()
{
	lm.OutputAscii();
}
