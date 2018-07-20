#include "fourthorder.h"

/*
Firstly, obviously the header file is included. 

This document contains the definition of the run function for the class that solves the membrane equation. 
The run function returns a double, which will be the energy of the surface. The first thing done is a surface refinemnet loop is 
written, with the KellyErrorEstimator. This checks to see if this is the first run through refinement to start, and if it is, then the 
surface Triangulation must be empty. Thus, the cell_mesh function is run to generate a coarse surface with the parameters we want. Then, 
the remove_anisotropy function is run on the surface with the golden ratio as its tolerance ofr size ratios between squares, and an allowance 
for number of iterations. This will slice up the more disproportionate cells of the coarse mesh into squarer pieces. Then, the system is solved for
w values, and the two-dimensional mesh is written into an eps file for viewing. Once the first run has been completed, the loop then estimates the error 
each cell of the program, and refines the sections with less error, while coarsening the sections with more error. Once this is done, the final 3d surface 
written to a .gpl file, and the energy is of the surface is calculated using the 
*/

double FourthOrder::run(double r1, double r2, double sep, double x, double y, double sigma, double kappa, double kappabar, int i){

	for(unsigned int refine_cycle = 0; refine_cycle < 6; ++refine_cycle){
		if(refine_cycle == 0){
			cell_mesh(r1, r2, sep, x, y, true);

			GridTools::remove_anisotropy(surface, 1.6180339887, 5);
		}

		else{
			Vector<float> estimated_error(surface.n_active_cells());
			KellyErrorEstimator<2>::estimate(doffer, QGauss<1>(3), typename FunctionMap<2>::type(), solution, estimated_error);

			GridRefinement::refine_and_coarsen_fixed_number(surface, estimated_error, .25, 0.05);
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