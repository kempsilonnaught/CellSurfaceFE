#include "TwoInclusionLift.h"

void SolveLaplacian::setup(){
	doffer.distribute_dofs(fe);
	DynamicSparsityPattern dynspar(doffer.n_dofs());
	DoFTools::make_sparsity_pattern(doffer, dynspar);
	sparsity_pattern.copy_from(dynspar);

	BIG_matrix.reinit(sparsity_pattern);

	solution.reinit(doffer.n_dofs());
	rhs.reinit(doffer.n_dofs());
}