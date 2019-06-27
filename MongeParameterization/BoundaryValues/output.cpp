//
// Created by Sam on 5/12/2019.
//
#include "forces_inclusions.h"


/*
This is the output function. Relatively simple, and solely for visualization purposes. This functions outputs a gnuplot file
that allows us to create a three dimensional image of the surface. This can either be viewed in gnuplot or be read into a mathematica file.
We have written a relatively basic mathematica program that parses gnuplot files, and from there writing a program to view the surface is trivial.
*/

void SimulateSurface::output(int i, int j){
    DataOut<2> data_out;

    data_out.attach_dof_handler(doffer);
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches();

    if(j = 1){
        std::ofstream out("gpls/positive/surface" + std::to_string(i) + "positive.gpl");
    	data_out.write_gnuplot(out);
    }else if(j = -1){
        std::ofstream out("gpls/negative/surface" + std::to_string(i) + "negative.gpl");
        data_out.write_gnuplot(out);
    }else{
        std::cout << "Error: Inclusion direction indicator (j) has an invalid value." << std::endl;
    }
