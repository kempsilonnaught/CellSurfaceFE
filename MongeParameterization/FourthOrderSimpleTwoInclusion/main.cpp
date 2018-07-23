#include "fourthorder.h"

/*
The main file is used to control the whole program. Thus, after including the header file,
the contructor and destructor are definedfor our FourthOrder class, which is defined in the header
file. The constructor initializes our fe object, which was declared in the header file,
with a numerical argument that indicates the order of the polynomials based on Gauss-Lobatto nodes.

Next is the main function. This function organizes the running of the program, and facilitates it 
being run in an iterative manner. Firstly, three arrays are declared. The first two are of type double, 
and will be where we write our data to. The Energy array will hold the Energy values for each run through the 
program, and the Separation array will hold the separation from each run through the program. The third array is 
an array of instantiations of the class for solving the problem. This is so when we loop to run the program multiple times
we don't have any data leftover from the last run due to the deconstructor not being called properly. The aforementioned was a 
large problem while writing the program. All three arrays are of an arbitrarily large size so long as they are large enough to hold all of 
the data for the given number of loops. A counting integer i is declared and set to 0, and a text file to record the energy and separation of
each run through the program is opened. 

The for loop allows us to incrementally move the inclusions away from each other on the bilayer, recording differences in energy. Ultimately,
the program overall produces the energy of the surface for a given separation of the inclusion. Thus, to take data, we loop over 
separations starting at 50, and going to 1000 at increments of 5. The radii of the inclusions are then declared and defined, as well as 
a few of the constants in our equation, such as surface tension(sigma), (kappa), and (kappabar).

*/

FourthOrder::FourthOrder() : fe(1){}
FourthOrder::~FourthOrder(){}

int main(){

	double Energy[1000];
	double Separation[1000];
	FourthOrder solve_instance[1000];
	unsigned int i = 0;

	std::ofstream energysep;
	energysep.open("energysep.txt");

	for(double sep = 50; sep <= 1000; sep += 25){
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
		std::cout << Energy[i] << std::endl;

		++i;
	}

	energysep.close();

	return 0;
}
