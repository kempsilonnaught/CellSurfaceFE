#include "fourthorder.h"

/*

*/

double FourthOrder::run(double r1, double r2, double sep, double x, double y, double sigma, double kappa, double kappabar, int i){

	for(unsigned int refine_cycle = 0; refine_cycle < 8; ++refine_cycle){
		if(refine_cycle == 0){
			cell_mesh(r1, r2, sep, x, y, true);

			GridTools::remove_anisotropy(surface, 1.6180339887, 8);
		}

		else{
			Vector<float> estimated_error(surface.n_active_cells());
			KellyErrorEstimator<2>::estimate(doffer, QGauss<1>(3), typename FunctionMap<2>::type(), solution, estimated_error);

			GridRefinement::refine_and_coarsen_fixed_number(surface, estimated_error, .30, 0.10);
			surface.execute_coarsening_and_refinement();
		}

		setup();
		assemble(sigma, kappa, kappabar);
		solve();


		std::ofstream out("twoDgrids/testgrid" + std::to_string(i) + ".eps");
		GridOut cell_mesho;
		cell_mesho.write_eps(surface, out);
	}

	std::cout << sep << std::endl;
	output(i);
	double energy = calcEnergy(sigma, kappa, kappabar);

	cell_mesh(r1, r2, sep, x, y, false);
	doffer.clear();
	surface.clear();

	return energy;
}