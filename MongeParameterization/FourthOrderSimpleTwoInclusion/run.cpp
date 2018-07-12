#include "fourthorder.h"

double FourthOrder::run(double r1, double r2, double sep, double x, double y){

	cell_mesh(r1, r2, sep, x, y, true);
	setup();
	assemble();
	solve();
	output();

	double energy = calcEnergy();

	cell_mesh(r1, r2, sep, x, y, false);

	return energy;
}