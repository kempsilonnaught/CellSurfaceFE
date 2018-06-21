#include "TwoInclusionLift.h"

/* int main(int argc, char* argv[]){
	if(argc != what ever I want){
		std::cerr << "";
		return 1; 
	}
	int num = atoi(argv[2]);
	double sep = atof(argv[3]);
}
*/

SolveLaplacian::SolveLaplacian(){}
SolveLaplacian::~SolveLaplacian(){}

int main(){
	SolveLaplacian mysys;
	mysys.cell_mesh();
}