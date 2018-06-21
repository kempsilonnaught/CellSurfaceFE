#ifndef TWO_INCLUSION_LIFT_SAMKEMP
#define TWO_INCLUSION_LIFT_SAMKEMP

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>


using namespace dealii;

class SolveLaplacian{
public:
	SolveLaplacian();
	~SolveLaplacian();

	void cell_mesh();
	void setup();
	void assemble();
	void solve();
	void output();

private:

	Triangulation<2> base_surface;
	Triangulation<2> hole_1;
	Triangulation<2> hole_2;
	Triangulation<2> remove_mesh;
	Triangulation<2> final_surface;
	/*DoFHandler<2> doffer;

	FESystem<2> fe;

	ConstraintMatrix constraint;

	SparsityPattern very_sparse;
	SparseMatrix<double> sparse_matrix; */

	Vector<double> later;
	Vector<double> latertwo;

};

#endif