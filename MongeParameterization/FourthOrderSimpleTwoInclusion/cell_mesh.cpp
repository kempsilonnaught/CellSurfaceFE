#include "fourthorder.h"

/*
Obviously the header file is included above. This file defines the cell_mesh function which
generates Triangulation in the surface object. The program takes as arguments the radii of the 
inclusions, the distance between them, and the overall x and y dimensions of the base manifold. 
This function also takes a boolean that indicates whether or not this is the first time through the cell_mesh 
function for this iteration of the program. Note that the run file runs this function twice. 

Firstly, the surface object is cleared of any remnant data, now that it has entered the scope where everything it references
is defined. Next two separate meshes are created, both being squares with a hole in the center, centered on the origin. The squares have side length
equal to the value selected for separation between the two inclusion, by defining half of a side length to be half of the separation. Then, one mesh
is moved in the positive x direction, and the other in the negative x direction, both by a value of half of the separation. Thus, the edges of each square 
meet at x = 0, and the two holes have a distance equal to the separation value between them. The two triangulations are then merged into one, and this 
triangulation is attached to the surface object. Then, two pointers to circular boundaries are created, with each circular boundary defined to be a circle of
radius equal to one of the inclusions, with the boundary centered on that inclusion. The boundaries are then set with an id. Next, we loop over all lines of all cells
in the surface, checking to see if they are at the boundary. If they are indeed on the boundary, we check to see if they are on the inclusion boundary, rather than the 
edge of the rectangular manifold. This is done by checking if the center of the line in question is less than or equal to a distance equal to the inclusion radius 
away from the center of the inclusion. If this is the case, the boundary id of the line is set to the circular boundary. This means that whenever we refine the mesh, these lines
divide in such a way that they become closer to being a circle. If the math regarding the center of the line and the radius does not make sense to you, I urge you to
draw a picture of the whole system, and draw an octogan and a circle either inscribed or circumscribed, preferably one of each, with the mesh stretching away from the nodes. 
Convince yourself that the center of the lines that are edges of the octogans must be less than or equal to the radius of the circle. Reread the code, and see that that is exactly 
what I have done here.

Nextly, 
*/

void FourthOrder::cell_mesh(double r1, double r2, double sep, double x, double y, bool first_run){
	surface.clear();

		Triangulation<2> inclusion_1;
		GridGenerator::hyper_cube_with_cylindrical_hole(inclusion_1, r1, sep/2, true);
		GridTools::shift(Point<2>(sep/2, 0), inclusion_1);  

		Triangulation<2> inclusion_2;
		GridGenerator::hyper_cube_with_cylindrical_hole(inclusion_2, r2, sep/2, true);
		GridTools::shift(Point<2>(-sep/2, 0), inclusion_2);

		GridGenerator::merge_triangulations(inclusion_1, inclusion_2, surface);

		HyperBallBoundary<2>* inclusion_boundary_1 = new HyperBallBoundary<2>(Point<2>(sep/2, 0), r1);
		surface.set_boundary(5, *inclusion_boundary_1); 

		for(unsigned int aa = 0; aa < 10; ++aa){
			typename Triangulation<2>::active_cell_iterator cell_1 = surface.begin_active(), endc_1 = surface.end();
				for(; cell_1 != endc_1; ++cell_1)
					for(unsigned int l_1 = 0; l_1 < GeometryInfo<2>::lines_per_cell; ++l_1){
						Point<2> edge_center_1 = cell_1 -> line(l_1) -> center(true, true);
							if(cell_1 -> line(l_1) -> at_boundary()){
								if(sqrt((std::pow((edge_center_1[0]-sep/2), 2))+std::pow(edge_center_1[1], 2)) <= (r1))
									cell_1 -> line(l_1) -> set_all_boundary_ids(5);
							}
					}
		}

		HyperBallBoundary<2>* inclusion_boundary_2 = new HyperBallBoundary<2>(Point<2>(-sep/2, 0), r2);
		surface.set_boundary(6, *inclusion_boundary_2); 
		
		for(unsigned int aa = 0; aa < 10; ++aa){
			typename Triangulation<2>::active_cell_iterator cell_2 = surface.begin_active(), endc_2 = surface.end();
				for(; cell_2 != endc_2; ++cell_2)
					for(unsigned int l_2 = 0; l_2 < GeometryInfo<2>::lines_per_cell; ++l_2){
						Point<2> edge_center_2 = cell_2 -> line(l_2) -> center(true, true);
							if(cell_2 -> line(l_2) -> at_boundary()){
								if(sqrt((std::pow((edge_center_2[0]+sep/2), 2))+std::pow(edge_center_2[1], 2)) <= (r2))
									cell_2 -> line(l_2) -> set_all_boundary_ids(6);
							}
					}
		}

		typename Triangulation<2>::active_cell_iterator 
		cell = surface.begin_active(), 
		endc = surface.end();
		for (; cell != endc; ++cell){
			for(unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i){
				Point<2> &v = cell -> vertex(i);
				if(v(0) == sep)
					v(0) = x/2; 
				if(v(0) == -sep)
					v(0) = -x/2;
				if(v(1) == sep/2)
					v(1) = y/2; 
				if(v(1) == -sep/2)
					v(1) = -y/2;              
			
			}

		}

		if(first_run == false){ 
			surface.clear();
			delete inclusion_boundary_1;
			delete inclusion_boundary_2;
		}

}

void FourthOrder::smoothness(){
	
}
		