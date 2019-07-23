#include "forces_inclusions.h"

void SimulateSurface::refine_mesh(double epsilon){
	Vector<float> estimated_error_per_cell(surface.n_active_cells());

	KellyErrorEstimator<2>::estimate(doffer, QGauss<1>(3), std::map<types::boundary_id, const Function<2> *>(), solution, estimated_error_per_cell);

	GridRefinement::refine_and_coarsen_fixed_number(surface, estimated_error_per_cell, 0.25, 0.03);

	surface.execute_coarsening_and_refinement();

	// for(unsigned int step=0; step<1; ++step){
	// 	Triangulation<2>::active_cell_iterator cell = surface.begin_active();
	// 	Triangulation<2>::active_cell_iterator endc = surface.end();
	// 	for(; cell!=endc; ++cell){
	// 		for(unsigned int l_1 = 0; l_1 < GeometryInfo<2>::lines_per_cell; ++l_1){
	// 			Point<2> edge_center_1 = cell -> line(l_1) -> center(true, true);
	// 			if(cell -> line(l_1) -> at_boundary()){
	// 				if(sqrt((std::pow((edge_center_1[0]-sep/2), 2))+std::pow(edge_center_1[1], 2)) <= (r1 + epsilon))
	// 					cell->set_refine_flag ();
	// 				break;
	// 			}
	// 		}
	// 		for(unsigned int l_2 = 0; l_2 < GeometryInfo<2>::lines_per_cell; ++l_2){
	// 			Point<2> edge_center_2 = cell -> line(l_2) -> center(true, true);
	// 				if(cell -> line(l_2) -> at_boundary()){
	// 					if(sqrt((std::pow((edge_center_2[0]+sep/2), 2))+std::pow(edge_center_2[1], 2)) <= (r2 +epsilon))
	// 						cell->set_refine_flag ();
	// 					break;
	// 				}
	// 		}
	// 	}
	// 	surface.execute_coarsening_and_refinement ();
	// }


}