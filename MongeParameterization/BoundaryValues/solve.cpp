//
// Created by Sam on 5/12/2019.
//
#include "forces_inclusions.h"

/*
Two functions are contained in this file. The solve function, which solves the matrix equation for w values, and the calcEnergy fucntion which calculates the total energy of the surface.

Firstly, the solve function which calculates our w values, which in turn give us our h values(technically, due to math trickery, and zeroes and ones, w = h).
We because the sparse matrix has all the rest of equation in it, the actual solve step is pretty simple. We set limits on the solver, giving it a maximum number of iterations
to reach convergence, and a level of error that should be considered close enough to count as adequate convergence. We then solve the equation, writing our w values to the solution vector.
Then as a control step for when we increase the solvers max number of iterations, we have the program print out how many steps it took to converge. This gives us a basis to check so we
don't accidentally set the lab computers on fire. Socks is squealing during this step as is.

Next, the calcEnergy function goes back and solves the orginal expression for Energy using the height values. Note that this is the final step in what is effectively a basic
calculus 1 minimzation technique. We had an expression for energy, took the derivative, set it equal to zero, solved for our variables at that, and now we are plugging our
variables back into our original Energy equation to get the minimum energy. We loop over all cells just like in the assemble function, only this time we loop using the Energy expression,
calculating the scalar value for the energy of every cell, and then adding them together. We store the total energy of the surface in a double and then return it.


*/

void SimulateSurface::solve(){
    SolverControl solver_control(1000000, 1e-11);
    SolverCG<> solver(solver_control);

    solver.solve(big_matrix, solution, rhs, PreconditionIdentity());

    std::cout << "   " << solver_control.last_step()
              << " CG iterations needed to obtain convergence."
              << std::endl;

}

double SimulateSurface::calcEnergy(double sigma, double kappa, double kappabar, double neumann_value){
    QGauss<2> quadrakilltwo(2);
    QGauss<1> face_quadrature_formula2(2);
    FEValues<2> fe_val(fe, quadrakilltwo, update_values | update_gradients | update_JxW_values | update_hessians);
    FEFaceValues<2> fe_face_val(fe, face_quadrature_formula2, update_values | update_quadrature_points | update_normal_vectors | update_hessians | update_JxW_values);


    unsigned int numdofs = fe.dofs_per_cell;
    unsigned int n_quadp = quadrakilltwo.size();
    unsigned int n_quadbound = face_quadrature_formula2.size();

    double energy = 0;
    double energy_surf = 0;
    double energy_bound = 0;

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
					energy_surf += ((kappa*((solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(trace(hess_i))*(trace(hess_j)))))*(fe_val.JxW(q))/2;
					energy_surf += ((kappabar*((hess_i[0][0])*(solution(local_dof_indices[i]))*(hess_j[1][1])*(solution(local_dof_indices[j])))))*(fe_val.JxW(q))/2;
                    energy_surf += ((kappabar*((hess_j[0][0])*(solution(local_dof_indices[j]))*(hess_i[1][1])*(solution(local_dof_indices[i])))))*(fe_val.JxW(q))/2;
					energy_surf += ((-kappabar*((solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(hess_i[0][1])*(hess_j[0][1]))))*(fe_val.JxW(q));
					energy_surf += ((sigma*((solution(local_dof_indices[i]))*(solution(local_dof_indices[j]))*(fe_val.shape_grad(i, q))*(fe_val.shape_grad(j, q)))))*(fe_val.JxW(q))/2;

				}
			}
		}


        for (unsigned int face_number = 0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
        	if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 5)){
            	fe_face_val.reinit (cell, face_number);
            	cell -> get_dof_indices(local_dof_indices);
        		for(unsigned int q = 0; q < n_quadbound; ++q){
                	for (unsigned int i=0; i<numdofs; ++i){
                		hess_i = fe_val.shape_hessian(i, q);
                    	energy_bound += -(solution(local_dof_indices[i]) * kappa * neumann_value * trace(hess_i) * fe_face_val.JxW(q));
                        energy_bound += -(solution(local_dof_indices[i]) * sigma * neumann_value * fe_face_val.JxW(q));
                    }
                }
            }

        for (unsigned int face_number = 0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
        	if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 6)){
            	fe_face_val.reinit (cell, face_number);
            	cell -> get_dof_indices(local_dof_indices);
        		for(unsigned int q = 0; q < n_quadbound; ++q){
                	for (unsigned int i=0; i<numdofs; ++i){
                		hess_i = fe_val.shape_hessian(i, q);
                    	energy_bound += -(solution(local_dof_indices[i]) * kappa * neumann_value * trace(hess_i) * fe_face_val.JxW(q));
                        energy_bound += -(solution(local_dof_indices[i]) * sigma * neumann_value * fe_face_val.JxW(q));
                    }
                }
            }

        energy += energy_bound + energy_surf;
    }

    return energy;
}
