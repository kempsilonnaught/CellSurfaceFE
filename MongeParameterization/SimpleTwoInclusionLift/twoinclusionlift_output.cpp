#include "TwoInclusionLift.h"

void SolveLaplacian::output(double s){
	DataOut<2> data_out;

	data_out.attach_dof_handler(doffer);
	data_out.add_data_vector(solution, "solution");

	data_out.build_patches();

	std::ofstream out("cell_forces" + std::to_string(s) + ".gpl");
	data_out.write_gnuplot(out);
}