#ifndef TWO_INCLUSION_LIFT_SAMKEMP
#define TWO_INCLUSION_LIFT_SAMKEMP

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>


#include <fstream>
#include <iostream>

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

class SolveLaplacian {

public :

	SolveLaplacian();
	~SolveLaplacian();

	void cell_mesh(double r1, double r2, double s, double x, double y, int n);
	void setup();
	void assemble();
	void solve();
	double calcEnergy();
	void output(double s);
	double run(double r1, double r2, double s, double x, double y);

private : 

	Triangulation<2> surface;
	double inclusion_rad_1;
	double inclusion_rad_2;
	double inclusion_separation;
	double surface_y;
	double surface_x; 
	
	DoFHandler<2> doffer;
	
	FE_Q<2> fe;
	
	SparsityPattern sparsity_pattern;
	
	SparseMatrix<double> big_matrix;
	
	Vector<double> solution;
	Vector<double> rhs;
	
};

#endif