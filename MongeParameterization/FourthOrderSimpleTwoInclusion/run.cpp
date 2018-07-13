#include "fourthorder.h"

double FourthOrder::run(double r1, double r2, double sep, double x, double y, double sigma, double kappa, double kappabar){

	cell_mesh(r1, r2, sep, x, y, true);

	std::ofstream out("twoDgrids/testgrid" + std::to_string(sep) + ".eps");
	GridOut cell_mesho;
	cell_mesho.write_eps(surface, out);

	setup();
	assemble(sigma, kappa, kappabar);
	solve();
	output();

	double energy = calcEnergy();

	cell_mesh(r1, r2, sep, x, y, false);

	doffer.clear();

	return energy;
}