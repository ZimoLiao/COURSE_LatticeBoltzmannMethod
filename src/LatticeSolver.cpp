#include "LatticeSolver.h"

void LatticeSolver::Ignore(int l)
{
	char buffer[101];
	for (int i = 0; i != l; i++) {
		fin.getline(buffer, 100);
	}
}

void LatticeSolver::InitParameter()
{
	/* broadcast parameters read from input file */
	int bufi[16] = { nx, ny, di, dj, \
		b[1], b[2] ,b[3] ,b[4] ,b[5] ,b[6] ,b[7] ,b[8] , \
		is_srt, is_mrt, \
		nlb,nle };
	double bufd[5] = { tau, omega_e, omega_ep, omega_q, omega_nu };

	MPI_Bcast(bufi, 16, MPI_INT, host, MPI_COMM_WORLD);
	MPI_Bcast(bufd, 5, MPI_DOUBLE, host, MPI_COMM_WORLD);

	if (rank != host) {
		nx = bufi[0];
		ny = bufi[1];
		di = bufi[2];
		dj = bufi[3];
		for (int i = 1; i != 9; i++) {
			b[i] = bufi[i + 3];
		}

		is_srt = bufi[12];
		is_mrt = bufi[13];

		nlb = bufi[14];
		nle = bufi[15];

		tau = bufd[0];
		omega_e = bufd[1];
		omega_ep = bufd[2];
		omega_q = bufd[3];
		omega_nu = bufd[4];
	}

	/* init local size */
	ni = nx / di;
	nj = ny / dj;
	ns = ni * rank;
	ne = ni * (rank + 1) - 1;

}

void LatticeSolver::InitBoundary()
{
	/* read all boundaries */
	fin.open("in/input");

	Ignore(25);
	int type, direct, is, ie, js, je;
	double rho, u, v;

	for (int ib = 0; ib != nlb; ib++) {
		fin >> type >> direct >> is >> ie >> js >> je >> rho >> u >> v;

		if (is == -1) { is = nx - 1; }
		if (js == -1) { js = ny - 1; }
		if (ie == -1) { ie = nx - 1; }
		if (je == -1) { je = ny - 1; }

		if (!(ie<ns || is>ne)) {
			is = max(0, is - ns);
			ie = min(ni - 1, ie - ns);

			LatticeBound lb(ib, type, direct, is, ie, js, je, rho, u, v);
			lp.InitBoundary(lb);
		}
	}

	fin.close();
	MPI_Barrier(MPI_COMM_WORLD);
}

void LatticeSolver::InitEntity()
{
	/* read entities */
	fin.open("in/input");

	Ignore(29 + nlb);
	double r, x, y, phi, ux0, uy0, uphi0;
	for (int ie = 0; ie != nle; ie++) {
		fin >> r >> x >> y >> phi >> ux0 >> uy0 >> uphi0;

		LatticeEntity newle(ie, r, x, y, phi, \
			double(rank)*ni, double(rank + 1.0)*ni - 1, 0, nj - 1, \
			ux0, uy0, uphi0);

		int exist = 0;
		if (newle.IsExist()) {
			exist = 1;
			lm.InitEntity(newle);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		int exist_host[64] = { 0 };
		MPI_Gather(&exist, 1, MPI_INT, exist_host, 1, MPI_INT, host, MPI_COMM_WORLD);

		if (rank == host) {
			exist = 0;
			for (int i = 0; i != size; i++) {
				exist += exist_host[i];
				// TODO: 给host中颗粒副本加入标识以便通信
			}
			le.push_back(newle);
		}
	}

	fin.close();
	MPI_Barrier(MPI_COMM_WORLD);
}

void LatticeSolver::StartOrder()
{
	int buf = 1;
	MPI_Status status;
	if (rank == host) {
		//MPI_Send(&buf, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);
	}
	else if (rank != size - 1) {
		MPI_Recv(&buf, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
		//MPI_Send(&buf, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);
	}
	else if (rank == size - 1) {
		MPI_Recv(&buf, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
	}
}

void LatticeSolver::EndOrder()
{
	int buf = 1;
	MPI_Status status;
	if (rank == host) {
		MPI_Send(&buf, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);
	}
	else if (rank != size - 1) {
		//MPI_Recv(&buf, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
		MPI_Send(&buf, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD);
	}
	else if (rank == size - 1) {
		//MPI_Recv(&buf, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
	}
}

LatticeSolver::LatticeSolver()
{
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* parameters initialization */
	if (rank == host) {
		fin.open("in/input");

		Ignore(5);
		fin >> nx >> ny;
		Ignore(3);
		fin >> di >> dj;
		Ignore(3);
		b[0] = 0;
		fin >> b[1] >> b[2] >> b[3] >> b[4] >> b[5] >> b[6] >> b[7] >> b[8];
		Ignore(5);
		fin >> is_srt >> tau;
		Ignore(2);
		fin >> is_mrt >> omega_e >> omega_ep >> omega_q >> omega_nu;
		Ignore(5);
		fin >> nlb;
		Ignore(4 + nlb);
		fin >> nle;

		fin.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	InitParameter();

	/* lattice units initialization */
	lm.Init(rank*ni, 0, ni, nj);

	// fluid-boundaries
	InitBoundary();

	// entities (particles)
	InitEntity();
}

void LatticeSolver::Calculate()
{
	lm.Force();
}

void LatticeSolver::PrintInfo()
{
	if (rank == host) {
		cout << "======================================\n";
		cout << "= LBSolver-Parallel                  =\n";
		cout << "======================================\n";
		cout << "= grid                               =\n";
		cout << "= ---------------------------------- =\n";
		cout << "=\tnx : " << setw(5) << nx << " |" << setw(5) << di << " | " << setw(5) << ni << "    =\n";
		cout << "=\tny : " << setw(5) << ny << " |" << setw(5) << dj << " | " << setw(5) << nj << "    =\n";
		cout << "= ---------------------------------- =\n";
		cout << "= number of entities :    " << setw(7) << nle << "    =\n";
		cout << "= number of boundaries :  " << setw(7) << nlb << "    =\n";
		cout << "======================================\n";
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void LatticeSolver::WriteUnit()
{
	if (rank == host) {
		/* calculate lattice type
		int* type = new int[nx*ny];
		for (int i = 0; i != nx * ny; i++) { type[i] = 0; }
		for (int ie = 0; ie != nle; ie++) {

		}
		delete[] type;
		*/

		/* write tecplot ascii format file */
		string fname = "out/unit_" + to_string(step) + ".dat";
		fout.open(fname);

		int nm = 0;
		for (int ie = 0; ie != nle; ie++) {
			nm += le[ie].get_nm();
		}

		// header
		fout << "TITLE     = Particles" << endl;
		fout << "FILETYPE  = FULL" << endl;
		fout << "VARIABLES = x y up vp" << endl;
		fout << "ZONE    F = point" << endl;
		fout << "        I = " << nm << endl;
		fout << "SOLUTIONTIME = " << step << endl;

		// particles
		for (int ie = 0; ie != nle; ie++) {
			for (int im = 0; im != le[ie].get_nm(); im++) {
				fout << le[ie].get_xm(im) << ' ' << le[ie].get_ym(im) << ' ' \
					<< le[ie].get_uxm(im) << ' ' << le[ie].get_uym(im) << endl;
			}
		}

		fout.close();
	}
}

void LatticeSolver::WriteFlow()
{
	string fname = "out/flow_" + to_string(rank) + "_" + to_string(step) + ".dat";
	fout.open(fname);

	lm.WriteAscii(fout, step);

	fout.close();
}

