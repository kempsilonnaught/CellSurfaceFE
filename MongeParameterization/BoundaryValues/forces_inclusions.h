//
// Created by Sam on 5/6/2019.
//

#ifndef BOUNDARYCONDITIONS_FORCES_INCLUSIONS_H
#define BOUNDARYCONDITIONS_FORCES_INCLUSIONS_H

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
#include <deal.II/base/table_handler.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>


#endif
/*
This is a header file to be included in all files that run this program. This holds the Class declaration for solving
the membrane shape equation, and holds all of the necessary base pieces for the program to run as it should. This file is
protected <describe the protecty thing>

Firstly, a class is created to solve any specific fourth order problem given
differing inputs. This way, the program can be run many times given differing positions
and sizes of inclusions, as well as differing surface size. All contained within a for loop
in the main function. While the use of a class may complicate things, it has proven overall
worth it for efficiency.

Below, a constructor and deconstructor are declared, as well as all of the functions necessary to
run a simulation. They are appropriately named. cell_mesh generates the surface, setup connects all
of the objects(dof_handler, surface, sparse_matrix....etc), assemble loops over all cells in surface
to create the matrices represented in the our theoretical function. (If you are reading this, you should
probably be familiar with the theory at this point, though not necessarily comfortable with it. If not, look up
the Helfrich Hamiltonian and the soap film/thin film problem, as well as some reading on finite element method).

Nextly, objects necessary throughout the program are declared, but not defined.
	- surface is declared as a two dimensional triangulation, this will hold the data describing the mesh of the bilayer.
	- r1 is a double holding the radius of the first inclusion.
	- r2 is a double holding the radius of the second inclusion.
	- sep is a double holding the distance between the inclusions.
	- x is a double indicating the length in of the base of the surface in the x direction.
	- y is a double indicating the length in of the base of the surface in the y direction.
	- doffer is a DoFHandler in two dimensions. This is used to handle and manipulate the degrees of freedom of each cell.
	- fe is our finite element object for our surface, which holds our finite element space.
	- big_matrix is the sparse matrix that will hold the solution to the delE equation
	- solution is a vector that will hold all of the w values calculated from setting delE equal to zero, thus finding the w values for a minimum energy.
	- rhs holds all of the right hand side values for the equation, which in our case, since we are minimizing, is a giant zero vector.

Above, all of the header files necessary to run the entire program are included, both from deal.II and the C++ libraries.
The default namespace is declared dealii, and references to the C libraries will use the classic "std::<command_being_used>" notation.

*/

using namespace dealii;

class SimulateSurface{

public :

    SimulateSurface();
    ~SimulateSurface();

    void cell_mesh(double r1, double r2, double sep, double x, double y, bool first_run);
    void setup();
    void assemble(double sigma, double kappa, double kappabar, double neumann_value_1, double neumann_value_2);
    void solve();
    double calcEnergy(double sigma, double kappa, double kappabar, double neumann_value_1, double neumann_value_2);
    void output(int i);
    double run(double r1, double r2, double sep, double x, double y, double sigma, double kappa, double kappabar, double neumann_value_1, double neumann_value_2, int i);

private :

    Triangulation<2> surface;
    double r1;
    double r2;
    double sep;
    double x;
    double y;

    DoFHandler<2> doffer;

    FE_Q<2> fe;

    SparsityPattern sparsity_pattern;

    SparseMatrix<double> big_matrix;

    Vector<double> solution;
    Vector<double> rhs;
};

