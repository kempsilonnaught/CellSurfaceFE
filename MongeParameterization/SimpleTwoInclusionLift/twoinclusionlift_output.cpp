#include "TwoInclusionLift.h"

void SolveLaplacian::output(){
	DataOut<2> data_out;

	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");

	data_out.build_patches();

	std::ofstream out("cell_forces.eps");
	data_out.write_eps(out);
}