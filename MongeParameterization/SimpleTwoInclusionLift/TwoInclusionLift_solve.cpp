#include "TwoInclusionLift.h"

void SolveLaplacian::solve(){
	SolverControl solver_control(1000, 1e-12);
	SolverCG<> solver(solver_control);

	solver.solve(BIG_matrix, solution, rhs, PreconditionIdentity());
}