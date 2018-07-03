#include "TwoInclusionLift.h"

void SolveLaplacian::solve(){

	SolverControl solver_control(100000, 1e-12);
	SolverCG<> solver(solver_control);

	solver.solve(big_matrix, solution, rhs, PreconditionIdentity());

//###################################################################

	QGauss<2> quadrakilltwo(2);
	FEValues<2> fe_val(fe, quadrakilltwo, update_gradients | update_JxW_values);

	const unsigned int numdofs = fe.dofs_per_cell;
	const unsigned int n_quadp = quadrakilltwo.size();

	double integrand;
	double delh_squared;

	std::vector<types::global_dof_index> local_dof_indices(numdofs);

	for(auto cell : doffer.active_cell_iterators()){
		fe_val.reinit(cell);
		for(unsigned int q = 0; q < n_quadp; ++q){
			delh_squared = 0;
			for(unsigned int i = 0; i < numdofs; ++i){
				for(unsigned int j = 0; j < numdofs; ++j){
					cell -> get_dof_indices(local_dof_indices); 
					for(unsigned int k = 0; k < numdofs; ++k)
						for(unsigned int l = 0; l < numdofs; ++l)
							delh_squared += (solution(local_dof_indices[k]))*(solution(local_dof_indices[l]))*(fe_val.shape_grad(i, q))*(fe_val.shape_grad(j, q));
				}
			}

			integrand += fe_val.JxW(q)*sqrt(1 + delh_squared);		
			
		}
	}

	std::cout << "ENERGY = " << integrand << std::endl;
}