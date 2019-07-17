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

    std::vector<bool> boundary_dofs5(doffer.n_dofs(), false);
    std::vector<bool> boundary_dofs6(doffer.n_dofs(), false);
    DoFTools::extract_boundary_dofs(doffer, ComponentMask(), boundary_dofs5, {5});
    DoFTools::extract_boundary_dofs(doffer, ComponentMask(), boundary_dofs6, {6});


    const unsigned int first_boundary_dof5 = std::distance(boundary_dofs5.begin(), std::find(boundary_dofs5.begin(), boundary_dofs5.end(), true));
    const unsigned int first_boundary_dof6 = std::distance(boundary_dofs6.begin(), std::find(boundary_dofs6.begin(), boundary_dofs6.end(), true));

    constraints5.clear();
    for(unsigned int i = first_boundary_dof5 + 1; i < doffer.n_dofs(); ++i){
        if(boundary_dofs5[i] == true){
            constraints5.add_line(i);
            constraints5.add_entry(i, first_boundary_dof5, 1);
        }
    }
    constraints5.close();

    constraints6.clear();
    for(unsigned int i = first_boundary_dof6 + 1; i < doffer.n_dofs(); ++i){
        if(boundary_dofs6[i] == true){
            constraints6.add_line(i);
            constraints6.add_entry(i, first_boundary_dof6, 1);
        }
    }
    constraints6.close();

    DynamicSparsityPattern dynspar(doffer.n_dofs());
    DoFTools::make_sparsity_pattern(doffer, dynspar);
    constraints5.condense(dynspar);
    constraints6.condense(dynspar);
    sparsity_pattern.copy_from(dynspar);

    big_matrix.reinit(sparsity_pattern);

    solution.reinit(doffer.n_dofs(), false);
    rhs.reinit(doffer.n_dofs(), false);

    std::cout << "Number of degrees of freedom: "
              << doffer.n_dofs()
              << std::endl;

}