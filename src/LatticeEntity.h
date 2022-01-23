#ifndef LATTICEENTITY_H_
#define LATTICEENTITY_H_

#include<math.h>

constexpr double pi = 3.1415926535897932;

/* finite-size particle (with IBM boundary markers)
*/
class LatticeEntity
{
	/* particle */
	double r;
	double x, y, theta;

	/* markers */
	int nm;
	double dtheta;
	double* xm, *ym;

	/* internal functions */
	void UpdateMarker();

public:
	LatticeEntity(double r, double x, double y);
	~LatticeEntity();

};

#endif // !LATTICEENTITY_H_