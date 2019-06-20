#include "forces_inclusions.h"

/*
This file defines the assemble function. The assemble function is the function that produces the expression for del E and populates our
sparse matrix so we can solve for our w vector and get our height values.

Firstly, it generates a Gaussian Quadrature of our system with 2 quadrature points in each space dimension. Then, we create an FEValues object
which holds the values that can be computed with the finite element method. It uses our fe object as well as our Gaussian Quadrature, and whenever it is 
initialized, it updates the values, gradients, JacobiansxW, and hessian matrices for whatever object it is called on. We then make objects that describe
the number of degrees of freedom and number of quadrature points for our cells. We create a small matrix to hold the values of an individual cell, as well 
as an object to serve the same purpose for the right side of our equation. However, we are doing a minimization problem, so our right side is just zero. Thus
all rhs matrices will be fully populated with zeroes. Nextly, we declare a numbering pattern for the degrees of freedom on a cell that differs from the 
numbering for the whole mesh overall. This means we can just loop over our nodes starting at 0 for every cell. Then, we declare 2 2-dimenional Tensors with
of size 2 in both dimensions. These will hold hessian values for our calculations. 

Nextly, we loop over all cells in our surface. The FE_Values object, called fe_time is reinitialized at the beginning of the loop for every cell. Thus, all
values, gradients, JxW's, and hessians are updated for each cell before calculation. The cell matrix and right hand side matrix are reset to 0 for the cell. 
Then we loop over all quadrature points and degrees of freedom, first calculating our hessians, and then adding each term of our del E equation to the cell size
matrix. Each term is multiplied by the JacobianxW to include the mapping between the real cell and reference cell. We now have the integral value for that cell. We then find the global
indices for that cell, and add all of the cells values to the sparse matrix for the whole mesh. We do the same for the right hand side, but again, the rhs is zero. 
This is done for every cell, resulting in a sparse matrix with the integral values for all dofs. We then set the boundary values such that the rectangular exterior boundary has a height of 0,
and the circular boundaries for the inclusions are at a height of 400. We then write those values to our sparse matrix. Finally, for debugging purposes, we output the number of non-zero matrix elements. 
*/

void SimulateSurface::assemble(double sigma, double kappa, double kappabar, double neumann_value){
	QGauss<2> quadrakill(2);
	QGauss<1> face_quadrature_formula(2);

	FEValues<2> fe_time(fe, quadrakill, update_values | update_gradients | update_JxW_values | update_hessians);
    FEFaceValues<2> fe_face_values(fe, face_quadrature_formula, update_values | update_quadrature_points | update_normal_vectors | update_hessians | update_JxW_values);

    const Solution<2> exactish_solution;

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrakill.size();
	const unsigned int n_face_q_points = face_quadrature_formula.size();
	
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
					lil_matrix(i, j) += (kappa*((trace(hess_i))*(trace(hess_j))*(fe_time.JxW(q_index)))/2);
					lil_matrix(i, j) += (kappabar*((hess_j[0][0]))*(hess_i[1][1])*(fe_time.JxW(q_index)))/2;
					lil_matrix(i, j) += (kappabar*((hess_i[0][0]))*(hess_j[1][1])*(fe_time.JxW(q_index)))/2;
					lil_matrix(i, j) += -(kappabar*((hess_i[0][1]))*(hess_j[0][1])*(fe_time.JxW(q_index)));
					lil_matrix(i, j) += (sigma*(((fe_time.shape_grad(i, q_index))*(fe_time.shape_grad(j, q_index))*(fe_time.JxW(q_index)))))/2;

					lil_rhs(i) += 0;
				}
			}
				
		}

		for (unsigned int face_number = 0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
        	if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 5)){
            	fe_face_values.reinit (cell, face_number);
            	for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point){
                	for (unsigned int i=0; i<dofs_per_cell; ++i){
                		hess_i = fe_face_values.shape_hessian(i, q_point);
                    	lil_rhs(i) += (neumann_value *
                                    fe_face_values.shape_value(i, q_point)*
                                    fe_face_values.JxW(q_point));
                    }
                }
            }

        for (unsigned int face_number = 0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
        	if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 6)){
            	fe_face_values.reinit (cell, face_number);
            	for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point){
                	for (unsigned int i=0; i<dofs_per_cell; ++i){
                		hess_i = fe_face_values.shape_hessian(i, q_point);
                    	lil_rhs(i) += (neumann_value *
                                    fe_face_values.shape_value(i, q_point)*
                                    fe_face_values.JxW(q_point));
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
	MatrixTools::apply_boundary_values(boundary_values, big_matrix, solution, rhs);

	std::cout << "Number of non-zero Sparse Matrix entries: " << big_matrix.n_nonzero_elements() << std::endl;
}
