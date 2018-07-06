#include "TwoInclusionLift.h"

double SolveLaplacian::run(double r1, double r2, double s, double x, double y){
		//std::cout << s << "run" << std::endl;
		cell_mesh(r1, r2, s, x, y, 0);
		setup();
		//std::cout << s << "setup" << std::endl;
		assemble();
		//std::cout << s << "assemble" << std::endl;
		solve();
		//force(sigma);
		output(s);
		//std::cout << "Blep" + std::to_string(s) << std::endl;

		double energy = calcEnergy();

		cell_mesh(r1, r2, s, x, y, 1);

		doffer.clear();
		surface.clear();

		//big_matrix.clear();
		return energy;
}


