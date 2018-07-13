#include "fourthorder.h"

FourthOrder::FourthOrder() : fe(1){}
FourthOrder::~FourthOrder(){}

int main(){

	double Energy[100];
	double Separation[100];
	FourthOrder solve_instance[1000];
	unsigned int i = 0;

	std::ofstream energysep;
	energysep.open("energysep.txt");

	//for(#; s <= #; s+= #){
		double r1 = 10;
		double r2 = 10;
		double x = 500;
		double y = 250;
		double sep = 100;
		double sigma = 1;
		double kappa = 1;
		double kappabar = 1;
		solve_instance[i].run(r1, r2, sep, x, y, sigma, kappa, kappabar);

		//Energy[i] = solveinstance[i].run();
		//Separation[i] = s;

		//energysep << Spearation[i] << " " << Energy[i] << std::endl;

		//++i
	//}

	energysep.close();

	return 0;
}