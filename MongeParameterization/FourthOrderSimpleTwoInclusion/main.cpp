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

	for(double sep = 50; sep <= 1000; sep += 5){
		double r1 = 10;
		double r2 = 10;
		double x = 4000;
		double y = 2000;
		double sigma = 1;
		double kappa = 1;
		double kappabar = 1;

		Energy[i] = solve_instance[i].run(r1, r2, sep, x, y, sigma, kappa, kappabar, i);
		Separation[i] = sep;
		energysep << Separation[i] << " " << Energy[i] << std::endl;

		++i;
	}

	energysep.close();

	return 0;
}