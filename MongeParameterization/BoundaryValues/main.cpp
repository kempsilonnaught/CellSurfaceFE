#include "forces_inclusions.h"

SimulateSurface::SimulateSurface() : fe(1){}
SimulateSurface::~SimulateSurface(){}

int main() {
    int max = 500;
    double* energy = new double[max];
    double* separation = new double[max];
    SimulateSurface* membrane = new SimulateSurface[max];

    unsigned int i = 0;
    std::ofstream energysep;
    energysep.open("energysep.txt");

    double radius_1 = 10;
    double radius_2 = 10;
    double sheet_x = 2000;
    double sheet_y = 1000;
    double sigma = 1;
    double kappa = 1;
    double kappabar = 1;

    for (double sep = 50; sep <= 1750; sep += 10) {
    	const double neumann_value = tan(3.14159265/4);
        energy[i] = membrane[i].run(radius_1, radius_2, sep, sheet_x, sheet_y, sigma, kappa, kappabar, neumann_value, i);
        separation[i] = sep;
        energysep << separation[i] << " " << energy[i] << std::endl;
        std::cout << energy[i] << std::endl;

        ++i;
    }

    //i = 0;
    //std::ofstream energysep2;
    //energysep2.open("energysep.txt");
//
    //for (double sep = 50; sep <= 750; sep += 10) {
    //	double neumann_value = -tan(3.14159265/4);
    //    energy[i] = membrane[i].run(radius_1, radius_2, sep, sheet_x, sheet_y, sigma, kappa, kappabar, neumann_value, i, "negative");
    //    separation[i] = sep;
    //    energysep2 << separation[i] << " " << energy[i] << std::endl;
    //    std::cout << energy[i] << std::endl;
//
    //    ++i;
    //}

    energysep.close();

    return 0;
}
