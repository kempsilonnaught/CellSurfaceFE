#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include "TwoInclusionLift.h"

void SolveLaplacian::setup(){
	doffer.initialize(surface, fe);

	doffer.distribute_dofs(fe);
	DynamicSparsityPattern dynspar(doffer.n_dofs());
	DoFTools::make_sparsity_pattern(doffer, dynspar);
	sparsity_pattern.copy_from(dynspar);

	big_matrix.reinit(sparsity_pattern);

	solution.reinit(doffer.n_dofs(), false);
	rhs.reinit(doffer.n_dofs(), false);
}