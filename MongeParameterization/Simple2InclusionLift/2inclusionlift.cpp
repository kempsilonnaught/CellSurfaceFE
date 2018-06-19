#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>

namespace TwoInclusionLift{
	using namespace dealii;


	template <int dim>
	class SolveLaplacian{
	public:
		SolveLaplacian();
		~SolveLaplacian();
		void runprog();

	private:
		void setup();
		void assemble();
		void solve();
		void output();

	Triangulation<dim> mesh;
	DoFHandler<dim> 
	}


}