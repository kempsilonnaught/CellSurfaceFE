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
	double r1 = 10;
	double r2 = 5;
	double s = 100;
	double x = 300;
	double y = 150;

	SolveLaplacian mysys;
	mysys.cell_mesh(r1, r2, s, x, y);
}