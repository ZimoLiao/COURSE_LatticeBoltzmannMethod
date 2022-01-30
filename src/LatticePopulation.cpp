#include "LatticePopulation.h"


void LatticePopulation::InitData(LatticeMoment& lm)
{
	FlushData();
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			CalculateEquilibrium(&data[Index(i, j)], &lm(i, j));
		}
	}

	// TODO velocity profile of inlet
	double Re = 100.0, D = 40.0, ubar, H = nj - 1.0, tau = 0.6;
	ubar = Re * (tau - 0.5) / 3.0 / D;
	for (int j = 0; j != nj; j++) {
		uin[j] = 6.0 * ubar * double(j) * (H - j) / H / H;
	}
}

void LatticePopulation::InitBoundary(LatticeBound& newlb)
{
	if (newlb.IsExist()) {
		lb.push_back(newlb);
		nlb++;
	}
}

void LatticePopulation::Stream()
{
	// 1 direction 
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 1] = data[9 * (sizej * (i - 1) + j) + 1];
		}
	}

	// 2 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 2] = data[9 * (sizej * (i)+(j - 1)) + 2];
		}
	}

	// 3 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 3] = data[9 * (sizej * (i + 1) + j) + 3];
		}
	}

	// 4 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 4] = data[9 * (sizej * (i)+(j + 1)) + 4];
		}
	}

	// 5 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 5] = data[9 * (sizej * (i - 1) + (j - 1)) + 5];
		}
	}

	// 6 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 6] = data[9 * (sizej * (i + 1) + (j - 1)) + 6];
		}
	}

	// 7 direction
	for (int i = 1; i <= ni; i++) {
		for (int j = 1; j <= nj; j++) {
			data[9 * (sizej * i + j) + 7] = data[9 * (sizej * (i + 1) + (j + 1)) + 7];
		}
	}

	// 8 direction
	for (int i = ni; i > 0; i--) {
		for (int j = nj; j > 0; j--) {
			data[9 * (sizej * i + j) + 8] = data[9 * (sizej * (i - 1) + (j + 1)) + 8];
		}
	}
}

void LatticePopulation::Boundary()
{
	for (int b = 0; b != nlb; b++) {
		switch (lb[b].GetType())
		{
		case 1: // NEBB-0
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {

					lb[b].CalculateNebb0(&data[Index(i, j)]);

					// corner: set reference density
					if (i == 0 || i == ni - 1) {
						if (j == 0) {
							double rho_ref = 0.0;
							int ind = Index(i, 1);
							for (int a = 0; a != 9; a++) {
								rho_ref += data[ind + a];
							}
							//if (i == 0) { rho_ref = 1.0; } // 设置入口参考密度为1.0
							lb[b].CalculateNebb0(&data[Index(i, j)], rho_ref);
						}
						else if (j == nj - 1) {
							double rho_ref = 0.0;
							int ind = Index(i, nj - 2);
							for (int a = 0; a != 9; a++) {
								rho_ref += data[ind + a];
							}
							//if (i == 0) { rho_ref = 1.0; } // 设置入口参考密度为1.0
							lb[b].CalculateNebb0(&data[Index(i, j)], rho_ref);
						}
					}
				}
			}
			break;
		case 2: // NEBB-V
			// TODO 偷懒直接写入口速度剖面了
			// TODO 这里将velocity profile中心处密度固定为参考密度【NO 该方法不行】
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					//if (j == nj / 2 || j == nj / 2 - 1) {
					//	lb[b].CalculateNebbV(&data[Index(i, j)], 1.0, uin[j], 0.0);
					//}
					//else {
					lb[b].CalculateNebbV(&data[Index(i, j)], uin[j], 0.0);
					//}
				}
			}
			break;
		case 3: // NEBB-P
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					lb[b].CalculateNebbP(&data[Index(i, j)]);
				}
			}
			break;
		case 4: // extrapolation (TODO 修改)
			int indin, indout;
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					indin = Index(i - 1, j);
					indout = Index(i, j);
					for (int a = 0; a != 9; a++) {
						data[indout + a] = data[indin + a];
					}
				}
			}
			break;

		case 5:
			// TODO 偷懒直接写入口速度剖面了
			for (int i = lb[b].GetIs(); i <= lb[b].GetIe(); i++) {
				for (int j = lb[b].GetJs(); j <= lb[b].GetJe(); j++) {
					lb[b].CalculateAbbV(&data[Index(i, j)], uin[j], 0.0);
				}
			}
			break;
		}
	}
}

void LatticePopulation::CollideSrt(LatticeMoment& lm)
{
	for (int i = 0; i != size; i++) {
		data[i] *= omegac;
	}
	double feq[9];
	int ind;
	for (int i = 0; i != ni; i++) {
		for (int j = 0; j != nj; j++) {
			CalculateEquilibrium(feq, &lm(i, j));
			ind = Index(i, j);
			for (int a = 0; a != 9; a++) {
				data[ind + a] += feq[a] * omega;
			}
		}
	}
}
