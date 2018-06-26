#include "TwoInclusionLift.h"

void SolveLaplacian::run(double r1, double r2, double s, double x, double y){
	cell_mesh(r1, r2, s, x, y);
	setup();
	assemble();
	solve();
	output();
}


