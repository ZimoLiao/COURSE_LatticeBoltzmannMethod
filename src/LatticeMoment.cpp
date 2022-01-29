#include "LatticeMoment.h"

void LatticeMoment::InitData()
{
	for (int ind = 0; ind != sizeij; ind++) {
		data[3 * ind] = 1.0;
		data[3 * ind + 1] = 0.0;
		data[3 * ind + 2] = 0.0;
	}
}

void LatticeMoment::InitDataShear()
{
	double x, y, ymid = (double(nj) - 1.0) / 2.0;

	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			x = x0 + i;
			y = y0 + j;

			data[Index(i, j, 0)] = 1.0;
			if (y > ymid) { data[Index(i, j, 1)] = 0.1; }
			else if (y < ymid) { data[Index(i, j, 1)] = -0.1; }
			else { data[Index(i, j, 1)] = 0.0; }

			if (abs(y - ymid) < (nj / 4.0 - 1.0)) {
				data[Index(i, j, 2)] = 0.01 * cos(4.0 * pi * x / (double(mpi_size) * ni)) \
					* cos(2 * pi * (y - ymid) / double(nj));
			}
			else {
				data[Index(i, j, 2)] = 0.0;
			}
		}
	}
}

void LatticeMoment::Update(LatticePopulation& lp)
{
	double u, v;
	diff = 0.0;

	int ind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			ind = Index(i, j);
			u = data[ind + 1];
			v = data[ind + 2];

			CalculateMoment(&data[ind], &lp(i, j));

			diff += sqrt((u - data[ind + 1]) * (u - data[ind + 1]) + \
				(v - data[ind + 2]) * (v - data[ind + 2]));
		}
	}
}

void LatticeMoment::WriteAscii(int step)
{
	std::ofstream fout;
	string fname = "out/flow_" + std::to_string(step) + \
		"_" + std::to_string(mpi_rank) + ".dat";
	fout.open(fname);

	// header
	fout << "FILETYPE  = FULL" << endl;
	fout << "VARIABLES = x, y, rho, u, v" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << ni << endl;
	fout << "        J = " << nj << endl;
	fout << "SOLUTIONTIME = " << step << endl;

	// flow variables (moments)
	int ind;
	for (int j = 0; j != nj; j++) {
		for (int i = 0; i != ni; i++) {
			ind = Index(i, j);

			fout << std::left << std::setw(8) << x0 + i \
				<< ' ' << std::left << std::setw(8) << y0 + j \
				<< ' ' << std::scientific << std::setprecision(12) << data[ind] \
				<< ' ' << std::scientific << std::setprecision(12) << data[ind + 1] \
				<< ' ' << std::scientific << std::setprecision(12) << data[ind + 2] \
				<< '\n';
		}
	}

	fout.close();
}
