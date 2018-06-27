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

SolveLaplacian::SolveLaplacian() : fe(1), doffer(surface){}
SolveLaplacian::~SolveLaplacian(){}

int main(){
	for(double s = 100; s <= 600; s+=25){

		double r1 = 15;
		double r2 = 15;
		double x = 1000;
		double y = 500;

		SolveLaplacian twoinclusionlift;

		std::cout << s << std::endl;

		twoinclusionlift.run(r1, r2, s, x, y);

	}

	return 0;
}	
