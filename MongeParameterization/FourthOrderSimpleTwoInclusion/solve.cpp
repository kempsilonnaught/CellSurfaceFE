#include "fourthorder.h"

/*

*/

void FourthOrder::solve(){
	SolverControl solver_control(50000000, 1e-15);
	SolverCG<> solver(solver_control);

	solver.solve(big_matrix, solution, rhs, PreconditionIdentity());

	std::cout << "   " << solver_control.last_step()
        << " CG iterations needed to obtain convergence."
        << std::endl;

}

double FourthOrder::calcEnergy(double sigma, double kappa, double kappabar){
	QGauss<2> quadrakilltwo(2);
	FEValues<2> fe_val(fe, quadrakilltwo, update_values | update_gradients | update_JxW_values | update_hessians);

	unsigned int numdofs = fe.dofs_per_cell;
	unsigned int n_quadp = quadrakilltwo.size();

	double energy = 0;

	std::vector<types::global_dof_index> local_dof_indices(numdofs);

	Tensor<2,2> hess_i;
	Tensor<2,2> hess_j;

	for(auto cell : doffer.active_cell_iterators()){
		fe_val.reinit(cell);
		cell -> get_dof_indices(local_dof_indices); 
		for(unsigned int q = 0; q < n_quadp; ++q){
			for(unsigned int i = 0; i < numdofs; ++i){
				for(unsigned int j = 0; j < numdofs; ++j){
					hess_i = fe_val.shape_hessian(i, q);
					hess_j = fe_val.shape_hessian(j, q);
					energy += ((kappa*((solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(trace(hess_i))*(trace(hess_j)))/2))*(fe_val.JxW(q));
					energy += ((kappabar*((hess_i[0][0])*(solution(local_dof_indices[i]))*(hess_j[1][1])*(solution(local_dof_indices[j])))))*(fe_val.JxW(q));
					energy += ((-kappabar*((solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(hess_i[0][1])*(hess_j[0][1]))))*(fe_val.JxW(q));
					energy += ((sigma*((solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(fe_val.shape_grad(i, q))*(fe_val.shape_grad(j, q)))/2))*(fe_val.JxW(q));
				}
			}
		}
	}

	return energy;
}