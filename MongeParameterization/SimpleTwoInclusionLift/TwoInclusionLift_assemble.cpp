
#include "TwoInclusionLift.h"

void SolveLaplacian::assemble(){
	QGauss<2> quadrakill(2);

	/*The third argument is represented by some behind-the-scenes C magic that effectively
	allows(for readability purposes) to say that the bar(|) can be read as "and". i.e. the third
	argument is update_values and update_gradients and update_JxW_values. Note that JxW represents 
	Jacobians crossed with quadrature weights.
	*/

	FEValues<2> fe_time(fe, quadrakill, update_values | update_gradients | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrakill.size();

	/* Note that because we are in 2d, dofs per vertex is one and a cell has 4
	vertices, so there are 4 degrees of freedom per cell
	*/ 
	FullMatrix<double> lil_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> lil_rhs(dofs_per_cell);

	std::vector<types::global_dof_index>local_dof_indices(dofs_per_cell);
	std::cout << "BLORPDOG" << std::endl;

	for(auto cell : doffer.active_cell_iterators()){
		fe_time.reinit(cell);
		lil_matrix = 0;
		lil_rhs = 0;
		for(unsigned int q_index = 0; q_index < n_q_points; ++q_index){ //This is looping over all phi i and j defined in our mesh for every cell in the surface
			for(unsigned int i = 0; i < dofs_per_cell; ++i)
				for(unsigned int j = 0; j < dofs_per_cell; ++j)
					lil_matrix(i, j) += ((fe_time.shape_grad(i, q_index))*(fe_time.shape_grad(j, q_index))*(fe_time.JxW(q_index)));
			for(unsigned int i = 0; i < dofs_per_cell; ++i)	
				lil_rhs(i) += (fe_time.shape_value(i, q_index)*1*fe_time.JxW(q_index));	
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
	MatrixTools::apply_boundary_values(boundary_values, big_matrix, solution, rhs);
}
