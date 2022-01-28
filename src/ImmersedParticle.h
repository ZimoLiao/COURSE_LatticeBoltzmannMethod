#ifndef IMMERSEDPARTICLE_H_
#define IMMERSEDPARTICLE_H_

#include<math.h>

#include"LatticeMoment.h"

class ImmersedParticle
{
	const double pi = 3.1415926535897932;

	/* geometry parameters */
	double r;
	int nm;
	double dphi;

	// domain limit
	double xmin, xmax, ymin, ymax;


	/* data */
	// variables in local coordinate
	double x, y, phi, ux, uy, uphi;
	bool exist;
	double* mx, *my, *mux, *muy;
	bool* mexist;

	// unforced velocity & force at markers
	double* mufx, *mufy, *mFx, *mFy;


	/* internal functions */
	inline bool IsExist(double mxi);

	double D(double  d);
	double Delta(int im, double i, double j);

public:
	friend class LatticeMoment;


	/* constructors & destructor */
	ImmersedParticle(\
		double x0, double y0, double ni, double nj, double r, \
		double gx, double gy, double gphi, \
		double ux, double uy, double uphi);
	~ImmersedParticle();


	/* functions */
	bool IsExist();

	void ForceMoment(LatticeMoment& lm);

	double Interpolation(int im, LatticeMoment& lm);
};


#endif // !IMMERSEDPARTICLE_H_