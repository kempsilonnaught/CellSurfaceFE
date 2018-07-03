#include <whatever the fuck is needed>

Vector<double> w_vector; // note to self - this exist - it is your solution vector


QGauss<2 or whatever> qaudkill(2 or whatever);
FEValues<2> fe_values(fe, quadkill, update_gradients | update_JxW_values);

//Something to get the global_dof from the local_dof

const unsigned int numdofs = fe.dofs_per_cell;
const unsigned int numquadpoints = quadkill.size();

double integrand = 0;
double delh_squared;

//Some dofhandler shit to iterate over cells
for(; cell!=endc; ++cell){
	fe_values.reinit(cell);
	for( unsinged int q = 0; q < numquadpoints; ++q){
		// lo op over the dof's in the cell
		delh_squared = 0;
		for( unsigned int i=0; i<numdofs, i++){
			for(unsigned int j=0; j<numdofs, j++){
				delh_squared += w_vector(global_dof(i))*w_vector(global_dof(j))*
				fe_values.shape_grad(i,q)*fe_values.shape_grad(j,q)
			}
		}
		// add to the integrand
		integrand += fe_values.JxW(q)*sqrt(1+delh_squared);
	}
}


