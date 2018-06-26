#ifndef TWO_INCLUSION_LIFT_SAMKEMP
#define TWO_INCLUSION_LIFT_SAMKEMP

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>


using namespace dealii;

/* 
Creating a Solve Laplacian class, for solving a Laplacian. Duh.
Public functions first have a constructor and deconstructor. 

System functions come next for the actual computation of the laplacian. 
cell_mesh creates the surface with inclusions, taking the radii of the inclusions,
	the distance between them, and the x and y dimensions of the surface.
setup <say what setup does>
assemble <say what assemble does>
solve <say what solve does>
output <say what output does, including file types and how to read them>

In the private class:
Triangulation<2> surface declares a two dimensional mesh called "surface"
5 Doubles are then declared to represent the inputs to cell_mesh.

*/

class SolveLaplacian{
public:
	SolveLaplacian();
	~SolveLaplacian();

	void cell_mesh(double r1, double r2, double s, double x, double y);
	void setup();
	void assemble();
	void solve();
	void output();

private:

	Triangulation<2> surface;
	double inclusion_rad_1;
	double inclusion_rad_2;
	double inclusion_separation;
	double surface_y;
	double surface_x; 

};

#endif