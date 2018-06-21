#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <iostream>
#include <fstream>
#include <cmath>

#include "TwoInclusionLift.h"


void SolveLaplacian::cell_mesh(){
	GridGenerator::hyper_rectangle(base_surface, Point<2>(-8, -4), Point<2>(8, 4), false);

	GridGenerator::hyper_ball(hole_1, Point<2>(-4, 0), double radius = 0.5);
	GridGenerator::hyper_ball(hole_2, Point<2>(4, 0), double radius = 0.5);
	GridGenerator::create_union_triangulation(hole_1, hole_2, remove_mesh);

	GridGenerator::create_triangulation_with_removed_cells(base_surface, remove_mesh, final_surface); 




	std::ofstream out("testgrid.eps");
	GridOut cell_mesh;
	cell_mesh.write_eps(final_surface, out);
} 