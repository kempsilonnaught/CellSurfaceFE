#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/base/logstream.h>

#include <iostream>
#include <fstream>
#include <cmath>

#include "TwoInclusionLift.h"

/* 
The first 5 lines pass the input values to the variable declared
	in the private section of the class.

*/

void SolveLaplacian::cell_mesh(double r1, double r2, double s, double x, double y, int n){
	//free(dealii::Triangulation<2>(surface));

	inclusion_rad_1 = r1;
	inclusion_rad_2 = r2;
	inclusion_separation = s;
	surface_y = y;
	surface_x = x;
	

	std::cout << inclusion_separation << "mesh.1" << std::endl;

	Triangulation<2> inclusion_1;
	GridGenerator::hyper_cube_with_cylindrical_hole(inclusion_1, inclusion_rad_1, inclusion_separation/2);
	GridTools::shift(Point<2>(inclusion_separation/2, 0), inclusion_1);	

	Triangulation<2> inclusion_2;
	GridGenerator::hyper_cube_with_cylindrical_hole(inclusion_2, inclusion_rad_2, inclusion_separation/2);
	GridTools::shift(Point<2>(-inclusion_separation/2, 0), inclusion_2);

	GridGenerator::merge_triangulations(inclusion_1, inclusion_2, surface);

	HyperBallBoundary<2>* inclusion_boundary_1 = new HyperBallBoundary<2>(Point<2>(inclusion_separation/2, 0), inclusion_rad_1);
	surface.set_boundary(1, *inclusion_boundary_1);	

	for(unsigned int aa = 0; aa < 10; ++aa){
		typename Triangulation<2>::active_cell_iterator cell_1 = surface.begin_active(), endc_1 = surface.end();
			for(; cell_1 != endc_1; ++cell_1)
				for(unsigned int l_1 = 0; l_1 < GeometryInfo<2>::lines_per_cell; ++l_1){
					Point<2> edge_center_1 = cell_1 -> line(l_1) -> center(true, true);
						if(cell_1 -> line(l_1) -> at_boundary()){
							if(sqrt((std::pow((edge_center_1[0]-inclusion_separation/2), 2))+std::pow(edge_center_1[1], 2)) <= (inclusion_rad_1))
								cell_1 -> line(l_1) -> set_all_boundary_ids(1);
						}
				}
	}


	//std::cout << s << "mesh.3" << std::endl;

	HyperBallBoundary<2>* inclusion_boundary_2 = new HyperBallBoundary<2>(Point<2>(-inclusion_separation/2, 0), inclusion_rad_2);
	surface.set_boundary(2, *inclusion_boundary_2);	
	

	
	for(unsigned int aa = 0; aa < 10; ++aa){
		typename Triangulation<2>::active_cell_iterator cell_2 = surface.begin_active(), endc_2 = surface.end();
			for(; cell_2 != endc_2; ++cell_2)
				for(unsigned int l_2 = 0; l_2 < GeometryInfo<2>::lines_per_cell; ++l_2){
					Point<2> edge_center_2 = cell_2 -> line(l_2) -> center(true, true);
						if(cell_2 -> line(l_2) -> at_boundary()){
							if(sqrt((std::pow((edge_center_2[0]+inclusion_separation/2), 2))+std::pow(edge_center_2[1], 2)) <= (inclusion_rad_2))
								cell_2 -> line(l_2) -> set_all_boundary_ids(2);
						}
				}
	}

	typename Triangulation<2>::active_cell_iterator 
	cell = surface.begin_active(), 
	endc = surface.end();
	for (; cell != endc; ++cell){
		for(unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i){
			Point<2> &v = cell -> vertex(i);
			if(-1e-5 < (v(1)-s/2) && (v(1)-s/2) < 1e-5)
				v(1) = y/2;
			if(-1e-5 < (v(1)+s/2) && (v(1)+s/2) < 1e-5)
				v(1) = -y/2;
			if(-1e-5 < (v(0)-s) && (v(0)-s) < 1e-5)
				v(0) = x/2;
			if(-1e-5 < (v(0)+s) && (v(0)+s) < 1e-5)
				v(0) = -x/2;
		}
	} 


	//std::cout << s << "mesh.5" << std::endl;

	surface.refine_global(3);

	std::ofstream out("testgrid.eps");
	GridOut cell_mesh;
	cell_mesh.write_eps(surface, out);
	std::cout << s << "mesh.2" << std::endl;

	GridTools::remove_anisotropy(surface, 1, 1);

	if(n == 1){	
		surface.clear();
		delete inclusion_boundary_1;
		delete inclusion_boundary_2;
	}

	std::cout << s << "mesh.3" << std::endl;

	//std::cout << s << "mesh.5" << std::endl;

	//std::cout << s << "mesh.6" << std::endl;

	/* std::ofstream out("testgrid.eps");
	GridOut cell_mesh;
	cell_mesh.write_eps(surface, out);
	*/

} 