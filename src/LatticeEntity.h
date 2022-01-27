#ifndef LATTICEENTITY_H_
#define LATTICEENTITY_H_

#include<mpi.h>
#include<math.h>
#include<iostream>

using std::cout;
using std::endl;

class LatticeEntity
{
	const double pi = 3.1415926535897932;

	/* geometry */
	double r;

	/* kinematics */
	double xg, yg; // global coordinates
	double xs, xe, ys, ye; // relative position
	double x, y, phi;
	double ux, uy, uphi;

	/* marker */
	int nm;
	double dphi;
	double* xm, *ym, *uxm, *uym, *uxf, *uyf;

	/* internal functions */
	double FuncD(double r);

public:

	/* index */
	int index_global;
	int num_rank;
	int rank[64] = { 0 };


	/* constructor & destructor */
	LatticeEntity(int index_global, double r, double xg, double yg, double phi, \
		double xs, double xe, double ys, double ye, \
		double ux0, double uy0, double uphi0);
	LatticeEntity(const LatticeEntity& le);
	~LatticeEntity();

	/* functions */
	bool IsExist();
	inline bool IsExist(int im);

	void CalculateU(int im, double i, double j, double u, double v);
	double CalculateFx(int im, double rho);
	double CalculateFy(int im, double rho);

	void Reset();

	int get_nm();
	double get_xm(int im);
	double get_ym(int im);
	double get_uxm(int im);
	double get_uym(int im);
};


#endif // !LATTICEENTITY_H_