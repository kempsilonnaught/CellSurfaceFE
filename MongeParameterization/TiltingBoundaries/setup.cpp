//
// Created by Sam on 5/12/2019.
//
#include "forces_inclusions.h"

/*
This defines the setup function. This is relatively straight forward. The dof_handler object,
named "doffer", is initialized with the surface that has already been generated in cell_mesh,
along with our finite element object. Then, a sparsity pattern is developed and applied to the
our matrix, such that we can ignore the vacancies in the matrix and decrease run time. Of course,
this is all done under the hood of the deal.II library. The solution vector and rhs vector are then
reinitialized with the proper number of components. Of course, because our rhs is zero, our rhs
vector will be entirely populated with zeroes in the assemble functino to come. Then, just for convenience
in debugging, we print out the number of degrees of freedom.
*/

void SimulateSurface::setup(){
    doffer.initialize(surface, fe);

    std::vector<bool> boundary_dofs(doffer.n_dofs(), false);
    DoFTools::extract_boundary_dofs(doffer, ComponentMask(), boundary_dofs);

    const unsigned int first_boundary_dof = std::distance(boundary_dofs.begin(), std::find(boundary_dofs.begin(), boundary_dofs.end(), true));

    constraints.clear();
    constraints.add_line(first_boundary_dof);
    for(unsigned int i = first_boundary_dof + 1; i < doffer.n_dofs(); ++i){
        if(boundary_dofs[i] == true){
            constraints.add_entry(first_boundary_dof, i, -1);
        }
    }
    constraints.close();

    DynamicSparsityPattern dynspar(doffer.n_dofs());
    DoFTools::make_sparsity_pattern(doffer, dynspar);
    sparsity_pattern.copy_from(dynspar);
    constraints.condense(dynspar);

    big_matrix.reinit(sparsity_pattern);

    solution.reinit(doffer.n_dofs(), false);
    rhs.reinit(doffer.n_dofs(), false);

    std::cout << "   Number of degrees of freedom: "
              << doffer.n_dofs()
              << std::endl;

}