#include "LatticeSolver.h"

void LatticeSolver::Ignore(int l)
{
	char buffer[100];
	for (int i = 0; i != l; i++) {
		fin.getline(buffer, 99);
	}
}

LatticeSolver::LatticeSolver()
{
	int err = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == host) {
		// read input
		fin.open("in/input");

		Ignore(2);
		fin >> case_name;
		Ignore(5);
		fin >> Nx >> Ny;
		Ignore(3);
		fin >> Dx >> Dy;
		Ignore(3);
		for (int i = 1; i != 9; i++) {
			fin >> B[i];
		}
		Ignore(5);
		fin >> SRT >> tau;
		Ignore(2);
		fin >> MRT >> omega_e >> omega_ep >> omega_q >> omega_nu;

		fin.close();

		// print header lines
		string coll;
		if (SRT) { coll = "SRT"; }
		else { coll = "MRT"; }
		cout << "===========================\n";
		cout << "= LBSolver-Parallel start =\n";
		cout << "=-------------------------=\n";
		cout << "= grid       " << setw(4) << Nx << setw(4) \
			<< "*" << setw(4) << Ny << " =" << endl;
		cout << "= partition  " << setw(4) << Dx << setw(4) \
			<< "*" << setw(4) << Dy << " =" << endl;
		cout << "= number of processor " << setw(3) << size << " =" << endl;
		cout << "=-------------------------=\n";
		cout << "= collision operator  " << coll << " =\n";
		cout << "===========================\n\n";

	}

	// broadcast parameters 
	int buffer1[14] = { Nx, Ny, Dx, Dy, B[1], B[2] ,B[3] ,B[4] ,B[5] ,B[6] ,B[7] ,B[8] ,SRT, MRT };
	double buffer2[5] = { tau, omega_e, omega_ep, omega_q, omega_nu };
	MPI_Bcast(buffer1, 14, MPI_INT, host, MPI_COMM_WORLD);
	MPI_Bcast(buffer2, 5, MPI_DOUBLE, host, MPI_COMM_WORLD);

	if (rank != host) {
		Nx = buffer1[0];
		Ny = buffer1[1];
		Dx = buffer1[2];
		Dy = buffer1[3];
		for (int i = 1; i != 9; i++) {
			B[i] = buffer1[i + 3];
		}
		tau = buffer2[0];
		omega_e = buffer2[1];
		omega_ep = buffer2[2];
		omega_q = buffer2[3];
		omega_nu = buffer2[4];
	}

	// input partition do not match the number of processors
	// TODO: 需要一个全局报错退出函数
	if (size != (Dx * Dy)) {
		if (rank == host) {
			cout << "ERROR | part_proc_mismatch\n";
		}
		return;
	}

	lm.InitParameter(tau, omega_e, omega_ep, omega_q, omega_nu);
	lp.InitParameter(tau, omega_e, omega_ep, omega_q, omega_nu);

	lp.InitConnection(rank, Dx, Dy, B);
}
