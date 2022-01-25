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

	/* index */
	int index_global;

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
	double* xm, *ym, *uxm, *uym;

	/* internal functions */
	inline bool IsExist(int im);

public:
	/* constructor & destructor */
	LatticeEntity(int index_global, double r, double xg, double yg, double phi, \
		double xs, double xe, double ys, double ye, \
		double ux0, double uy0, double uphi0);
	LatticeEntity(const LatticeEntity& le);
	~LatticeEntity();

	/* functions */
	bool IsExist();
};


#endif // !LATTICEENTITY_H_