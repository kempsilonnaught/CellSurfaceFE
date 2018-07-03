#include <cstdlib>
#include <cmath>

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
SolveLaplacian::SolveLaplacian() : doffer(surface), fe(1){}
SolveLaplacian::~SolveLaplacian(){}

int main(){

	/* Takes Command Line Arguments:
	*/

	//int argc, char* argv[]
	//if(argc != 2){
	//	std::cout << "WRONG!" << std::endl;
	//	return 1;
	//}
	//double s = atof(argv[1]);
	
	/* Single Input: 
	*/

	//double sigma = 1;
	double r1 = 15;
	double r2 = 15;
	double x = 1000;
	double y = 500;
	double s = 300;
	
	SolveLaplacian twoinclusionlift;
	twoinclusionlift.run(r1, r2, s, x, y);

	/* Loops over increasing s:
	

	SolveLaplacian twoinclusionlift[1000];
	unsigned int i = 0;

	for(double s = 100; s < 750; s+=25){
		//double sigma = 1;
		double r1 = 15;
		double r2 = 15;
		double x = 1000;
		double y = 500;

		std::cout << s << std::endl;

		lift_i = twoinclusionlift[i];

		lift_i.run(r1, r2, s, x, y);
		++i;
	}
	*/

	return 0;
}	
