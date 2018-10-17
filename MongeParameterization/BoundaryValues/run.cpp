#include "fourthorder.h"

/*
Firstly, obviously the header file is included. 

This document contains the definition of the run function for the class that solves the membrane equation. Ultimnately, this file simply contains the instructions
to calculate our energy. The method is just like a simple calculus 1 extremum problem. You take the derivative of an expression, set it equal to zero, solve for your variable,
and then use that variable to solve your original equation for the extreme value. That is exactly what we have done here, just with a slightly more complicated expression, 
and thickly disguised with C++ and deal.II jargon. 

Firstly, we have to create a surface and setup our fe tools. We call the cell mesh function with specified inclusion radii, a rectangle size in terms of x and y, the separation between inclusions, and a boolean 
indicating whether or not this is the first time we have called it. This generates a mesh for the surface with specified size, inclusion requirements, and separation. We then refine all of the cells twice, and then 
refine just the cells around the inclusions another single time. Then we remove anistropy and hanging nodes. This does not make the surface uniform, but does insure that the size and shape of the cells change smoothly so 
we do not have abrupt changes in values. We the run the setup function which applies the finite element method to our surface. 

Next, we run the assemble function. We already have an expression for Energy of which we have taken the derivative which we call del E. We applied the small gradient approximation, and threw away some negligibly small
terms. These terms will come back when we add more inclusions to our system, as they will cease to be negligible at that point. The assemble function applies our equation for del E to the surface with the finite element method. 
From the assemble function we get a matrix of values that describes our equation thoroughly enough that we can solve for w with this matrix and the rhs values. 

This occurs in the solve fucntion, which we call next. This function sets the right hand side equal to zero, and sets the rest of the expression equal to the rhs. Then,
it solves for a vector of w values using the matrix found in the assemble function. For bookkeeping purposes, we output some useful information about the cells in the mesh. Then we output the two dimensional mesh 
as an eps as well as the separation. Then we run the output function which makes a gnuplot file of the three dimensional final surface. 

Finally, we run the calcEnergy function, which plugs our w values back into the original expression for Energy and solves for the minimum energy. This is stored in a double, and returned by this function at the end for the 
main function to store in our data. Finally, the cell_mesh is run again with the boolean equal to false, and it will delete all pointers and references it contains, and the rest of our objects will be cleared of data
so that they are ready for the next instantiation of the class. Typically this should be done in the deconstructor automatically when we go out of scope of the class instance, but deal.II's use of constants and statics made 
this complicated. 

Thus, this function has run the necessary functions to get energy.
*/

double FourthOrder::run(double r1, double r2, double sep, double x, double y, double sigma, double kappa, double kappabar, int i){

	cell_mesh(r1, r2, sep, x, y, true);
	surface.refine_global(2);

	for(unsigned int step=0; step<1; ++step){
		Triangulation<2>::active_cell_iterator cell = surface.begin_active();
		Triangulation<2>::active_cell_iterator endc = surface.end();
		for(; cell!=endc; ++cell){
			for(unsigned int l_1 = 0; l_1 < GeometryInfo<2>::lines_per_cell; ++l_1){
				Point<2> edge_center_1 = cell -> line(l_1) -> center(true, true);
				if(cell -> line(l_1) -> at_boundary()){
					if(sqrt((std::pow((edge_center_1[0]-sep/2), 2))+std::pow(edge_center_1[1], 2)) <= (r1))
						cell->set_refine_flag ();
						break;
				}
			}
			for(unsigned int l_2 = 0; l_2 < GeometryInfo<2>::lines_per_cell; ++l_2){
				Point<2> edge_center_2 = cell -> line(l_2) -> center(true, true);
					if(cell -> line(l_2) -> at_boundary()){
						if(sqrt((std::pow((edge_center_2[0]+sep/2), 2))+std::pow(edge_center_2[1], 2)) <= (r2))
							cell->set_refine_flag ();
							break;
					}
			}
		}
		surface.execute_coarsening_and_refinement ();
	}
	GridTools::remove_anisotropy(surface, 1.6180339887, 2);
	GridTools::remove_hanging_nodes(surface, false, 20);

	setup();
	assemble(sigma, kappa, kappabar);
	solve();

	std::cout << "   Number of active cells: "
	<< surface.n_active_cells()
	<< std::endl
	<< "   Total number of cells: "
	<< surface.n_cells()
	<< std::endl;

	std::ofstream out("twoDgrids/testgrid" + std::to_string(i) + ".eps");
	GridOut cell_mesho;
	cell_mesho.write_eps(surface, out);

	std::cout << sep << std::endl;
	output(i);
	double energy = calcEnergy(sigma, kappa, kappabar);

	cell_mesh(r1, r2, sep, x, y, false);
	doffer.clear();
	surface.clear();

	return energy;
}
