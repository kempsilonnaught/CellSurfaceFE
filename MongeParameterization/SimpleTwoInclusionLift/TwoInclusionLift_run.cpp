#include "TwoInclusionLift.h"

void SolveLaplacian::run(double r1, double r2, double s, double x, double y){
		//std::cout << s << "run" << std::endl;
		cell_mesh(r1, r2, s, x, y, 0);
		std::cout << s << "mesh" << std::endl;
		setup();
		//std::cout << s << "setup" << std::endl;
		assemble();
		//std::cout << s << "assemble" << std::endl;
		solve();
		//force(sigma);
		output();
		//std::cout << "Blep" + std::to_string(s) << std::endl;
		cell_mesh(r1, r2, s, x, y, 1);
		
		doffer.clear();
		surface.clear();
		//big_matrix.clear();
}


