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

SolveLaplacian::SolveLaplacian() : fe(1), dof_handler(surface){}
SolveLaplacian::~SolveLaplacian(){}

int main(){
	double r1 = 25;
	double r2 = 15;
	double s = 600;
	double x = 1000;
	double y = 500;

	SolveLaplacian twoinclusionlift;
	twoinclusionlift.run();

	return 0;
}