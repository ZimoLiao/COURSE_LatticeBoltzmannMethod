#include "LatticeSolver.h"

void LatticeSolver::Ignore(ifstream & fin, int l)
{
	char buffer[101];
	for (int i = 0; i != l; i++) {
		fin.getline(buffer, 100);
	}
}

LatticeSolver::LatticeSolver()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* read configuration */
	ifstream fin;
	fin.open("in/flow.init"); // TODO: 改成主节点读入再广播

	int di, dj;

	Ignore(fin, 4);

	fin >> nstep >> nwrite;
	Ignore(fin, 5);

	fin >> nx >> ny;
	Ignore(fin, 3);

	fin >> di >> dj;
	ni = nx / di;
	nj = ny / dj;
	Ignore(fin, 3);

	fin >> tg[1] >> tg[2] >> tg[3] >> tg[4];
	if (rank != 0) { tg[3] = 1; }
	if (rank != size - 1) { tg[1] = 1; }
	Ignore(fin, 5);

	fin >> is_srt >> tau; // TODO: MRT实现尚未部署
	Ignore(fin, 7);

	fin >> nlb_total;

	fin.close();


	/* moment & population initialization */
	lm.InitGeom(rank*ni, 0, ni, nj, tg);
	lm.InitData();

	lp.InitGeom(rank*ni, 0, ni, nj, tg);
	lp.InitData(lm);
	lp.InitParam(tau);


	/* read fluid-boundaries */
	int type, direct, is, ie, js, je;
	double rho, u, v;

	fin.open("in/flow.init");

	Ignore(fin, 29);
	for (int i = 0; i != nlb_total; i++) {
		fin >> type >> direct >> is >> ie >> js >> je >> rho >> u >> v;
		Ignore(fin, 1);

		if (is < 0) { is += nx; }
		if (ie < 0) { ie += nx; }
		if (js < 0) { js += ny; }
		if (je < 0) { je += ny; }

		LatticeBound newlb(rank, ni, nj, type, direct, is, ie, js, je, rho, u, v);
		lp.InitBoundary(newlb);
	}

	fin.close();
}
