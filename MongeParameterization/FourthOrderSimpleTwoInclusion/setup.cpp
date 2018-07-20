#include "fourthorder.h"

/*

*/

void FourthOrder::setup(){
	doffer.initialize(surface, fe);

	DynamicSparsityPattern dynspar(doffer.n_dofs());
	DoFTools::make_sparsity_pattern(doffer, dynspar);
	sparsity_pattern.copy_from(dynspar);

	big_matrix.reinit(sparsity_pattern);

	solution.reinit(doffer.n_dofs(), false);
	rhs.reinit(doffer.n_dofs(), false);


	std::cout << "   Number of degrees of freedom: "
        << doffer.n_dofs()
        << std::endl;

}