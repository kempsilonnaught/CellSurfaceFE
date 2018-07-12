#include <deal.II/grid/grid_refinement.h>
#include "TwoInclusionLift.h"

double SolveLaplacian::run(double r1, double r2, double s, double x, double y){

	for(unsigned int refine_cycle = 0; refine_cycle < 9; ++refine_cycle){
		//std::cout << s << "run" << std::endl;
		if(refine_cycle == 0){
			cell_mesh(r1, r2, s, x, y, 0);

			surface.refine_global(2);

		}

		else{
			Vector<float> estimated_error(surface.n_active_cells());
			KellyErrorEstimator<2>::estimate(doffer, QGauss<1>(3), typename FunctionMap<2>::type(), solution, estimated_error);

			GridRefinement::refine_and_coarsen_fixed_number(surface, estimated_error, .25, 0.05);
			surface.execute_coarsening_and_refinement();

		}
	

		setup();
		//std::cout << s << "setup" << std::endl;
		assemble();
		//std::cout << s << "assemble" << std::endl;
		solve();

		std::ofstream out("twoDgrids/testgrid" + std::to_string(s) + ".eps");
		GridOut cell_mesho;
		cell_mesho.write_eps(surface, out);

		std::cout << refine_cycle << std::endl;
	}


			
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


