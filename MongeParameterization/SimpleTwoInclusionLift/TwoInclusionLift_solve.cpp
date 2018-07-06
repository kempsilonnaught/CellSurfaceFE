#include "TwoInclusionLift.h"

void SolveLaplacian::solve(){

	SolverControl solver_control(100000, 1e-12);
	SolverCG<> solver(solver_control);

	solver.solve(big_matrix, solution, rhs, PreconditionIdentity());
}

/*###################################################################

	QGauss<2> quadrakilltwo(2);
	FEValues<2> fe_val(fe, quadrakilltwo, update_gradients | update_JxW_values);

	unsigned int numdofs = fe.dofs_per_cell;
	unsigned int n_quadp = quadrakilltwo.size();

	double energy = 0;
	double delh_squared;

	std::vector<types::global_dof_index> local_dof_indices(numdofs);

	for(auto cell : doffer.active_cell_iterators()){
		fe_val.reinit(cell);
		cell -> get_dof_indices(local_dof_indices); 
		for(unsigned int q = 0; q < n_quadp; ++q){
			delh_squared = 0;
			for(unsigned int i = 0; i < numdofs; ++i){
				for(unsigned int j = 0; j < numdofs; ++j){
					delh_squared += (solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(fe_val.shape_grad(i, q))*(fe_val.shape_grad(j, q));
				}
			}

			energy += fe_val.JxW(q)*(sqrt(1 + delh_squared));		
			
		}
	}

	std::cout << "ENERGY = " << energy << std::endl;
}
*/

double SolveLaplacian::calcEnergy() {
	QGauss<2> quadrakilltwo(2);
	FEValues<2> fe_val(fe, quadrakilltwo, update_gradients | update_JxW_values);

	unsigned int numdofs = fe.dofs_per_cell;
	unsigned int n_quadp = quadrakilltwo.size();

	double energy = 0;
	double delh_squared;

	std::vector<types::global_dof_index> local_dof_indices(numdofs);

	for(auto cell : doffer.active_cell_iterators()){
		fe_val.reinit(cell);
		cell -> get_dof_indices(local_dof_indices); 
		for(unsigned int q = 0; q < n_quadp; ++q){
			delh_squared = 0;
			for(unsigned int i = 0; i < numdofs; ++i){
				for(unsigned int j = 0; j < numdofs; ++j){
					delh_squared += (solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(fe_val.shape_grad(i, q))*(fe_val.shape_grad(j, q));
				}
			}

			energy += fe_val.JxW(q)*(sqrt(1 + delh_squared));		
			
		}
	}

	
	return energy;
}
