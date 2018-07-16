#include "fourthorder.h"

/*

*/

void FourthOrder::output(int i){
	DataOut<2> data_out;

	data_out.attach_dof_handler(doffer);
	data_out.add_data_vector(solution, "solution");

	data_out.build_patches();

	std::ofstream out("gpls/surface" + std::to_string(i) + ".gpl");
	data_out.write_gnuplot(out); 
}