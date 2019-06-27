#include "forces_inclusions.h"

SimulateSurface::SimulateSurface() : fe(2){}
SimulateSurface::~SimulateSurface(){}

int main() {
    int max = 500;
    double* energy = new double[max];
    double* separation = new double[max];
    SimulateSurface* membrane = new SimulateSurface[max];

    unsigned int i = 0;
    int j;
    std::ofstream energysep;
    energysep.open("energysep.txt");

    double radius_1 = 10;
    double radius_2 = 10;
    double sheet_x = 2250;
    double sheet_y = 1250;
    double sigma = 1;
    double kappa = 1;
    double kappabar = 1;

    for (double sep = 50; sep <= 1500; sep += 5) {
        j = 1;
    	const double neumann_value_1 = tan(3.14159265/4);
        const double neumann_value_2 = tan(3.14159265/4);
        energy[i] = membrane[i].run(radius_1, radius_2, sep, sheet_x, sheet_y, sigma, kappa, kappabar, neumann_value_1, neumann_value_2, i, j);
        separation[i] = sep;
        energysep << separation[i] << " " << energy[i] << std::endl;
        std::cout << energy[i] << std::endl;

        ++i;
    }

    energysep.close();

    i = 0;
    std::ofstream energysep2;
    energysep2.open("energysep2.txt");

    for (double sep2 = 50; sep <= 1500; sep += 5) {
        j = -1;
        const double neumann_value_3 = tan(3.14159265/4);
        const double neumann_value_4 = -tan(3.14159265/4);
        energy[i] = membrane[i].run(radius_1, radius_2, sep2, sheet_x, sheet_y, sigma, kappa, kappabar, neumann_value_3, neumann_value_4, i, j);
        separation[i] = sep2;
        energysep2 << separation[i] << " " << energy[i] << std::endl;
        std::cout << energy[i] << std::endl;

        ++i;
    }

    energysep2.close();

    return 0;
}
