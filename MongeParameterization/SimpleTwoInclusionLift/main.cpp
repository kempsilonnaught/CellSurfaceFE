#include <cstdlib>
#include <cmath>
#include <fstream>

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


SolveLaplacian::SolveLaplacian() : fe(1){}
SolveLaplacian::~SolveLaplacian(){
	
}

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
	//double r1 = 15;
	//double r2 = 15;
	//double x = 1000;
	//double s = 300;
	
	//SolveLaplacian twoinclusionlift;

	//twoinclusionlift.run(r1, r2, s, x, y);

	/* Loops over increasing s:
	*/

	double Energy[50];
	double Separation[50];

	SolveLaplacian twoinclusionlift[1000];
	unsigned int i = 0;

	std::ofstream energysep;
	energysep.open("energysep.txt");
	
	for(double s = 150; s <= 900; s+=15){
		//double sigma = 1;
		double r1 = 15;
		double r2 = 15;
		double x = 1000;
		double y = 500;

		Energy[i] = twoinclusionlift[i].run(r1, r2, s, x, y);
		Separation[i] = s;
		std::cout << Energy[i] << std::endl;
		std::cout << Separation[i] << std::endl;
		
		energysep << Separation[i] << " " << Energy[i] << std::endl;

		++i;
	}
	energysep.close();

	return 0;
}	
