#include "TwoInclusionLift.h"

void SolveLaplacian::solve(){
	SolverControl solver_control(100000, 1e-12);
	SolverCG<> solver(solver_control);

	solver.solve(big_matrix, solution, rhs, PreconditionIdentity());
}