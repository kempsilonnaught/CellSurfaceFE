#include "fourthorder.h"

/*

*/

void FourthOrder::assemble(double sigma, double kappa, double kappabar){
	QGauss<2> quadrakill(2);

	FEValues<2> fe_time(fe, quadrakill, update_values | update_gradients | update_JxW_values | update_hessians);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrakill.size();
	
	FullMatrix<double> lil_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> lil_rhs(dofs_per_cell);
	std::vector<types::global_dof_index>local_dof_indices(dofs_per_cell);

	Tensor<2,2> hess_i;
	Tensor<2,2> hess_j;

	for(auto cell : doffer.active_cell_iterators()){
		fe_time.reinit(cell);
		lil_matrix = 0;
		lil_rhs = 0;
		for(unsigned int q_index = 0; q_index < n_q_points; ++q_index){ //This is our big system integral
			for(unsigned int i = 0; i < dofs_per_cell; ++i){
				for(unsigned int j = 0; j < dofs_per_cell; ++j){
					hess_i = fe_time.shape_hessian(i, q_index);
					hess_j = fe_time.shape_hessian(j, q_index);
					lil_matrix(i, j) += (kappa*((trace(hess_i))*(trace(hess_j))*(fe_time.JxW(q_index))/2));
					lil_matrix(i, j) += (kappabar*((hess_j[0][0]))*(hess_i[1][1])*(fe_time.JxW(q_index)));
					lil_matrix(i, j) += (kappabar*((hess_i[0][0]))*(hess_j[1][1])*(fe_time.JxW(q_index)));
					lil_matrix(i, j) += (kappabar*(-(hess_i[0][1]))*(hess_j[0][1])*(fe_time.JxW(q_index)));
					lil_matrix(i, j) += (sigma*(((fe_time.shape_grad(i, q_index))*(fe_time.shape_grad(j, q_index))*(fe_time.JxW(q_index)))/2));
				}
			}
				
		}

		cell -> get_dof_indices(local_dof_indices);
		for(unsigned int i = 0; i < dofs_per_cell; ++i)
			for(unsigned int j = 0; j < dofs_per_cell; ++j)
				big_matrix.add(local_dof_indices[i], local_dof_indices[j], lil_matrix(i, j));
		for(unsigned int i = 0; i < dofs_per_cell; ++i)
			rhs(local_dof_indices[i]) += lil_rhs(i);
	}

	std::map<types::global_dof_index, double> boundary_values;
	VectorTools::interpolate_boundary_values(doffer, 0, ZeroFunction<2>(), boundary_values);
	VectorTools::interpolate_boundary_values(doffer, 5, ConstantFunction<2>(100), boundary_values);
	VectorTools::interpolate_boundary_values(doffer, 6, ConstantFunction<2>(100), boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, big_matrix, solution, rhs);

	std::cout << "Number of non-zero Sparse Matrix entries: " << big_matrix.n_nonzero_elements() << std::endl;

}
