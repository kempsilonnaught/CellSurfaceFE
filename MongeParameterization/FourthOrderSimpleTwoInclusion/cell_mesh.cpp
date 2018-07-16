#include "fourthorder.h"

/*

*/

void FourthOrder::cell_mesh(double r1, double r2, double sep, double x, double y, bool first_run){
	surface.clear();

    Triangulation<2> inclusion_1;
    GridGenerator::hyper_cube_with_cylindrical_hole(inclusion_1, r1, sep/2, true);
    GridTools::shift(Point<2>(sep/2, 0), inclusion_1);  

    Triangulation<2> inclusion_2;
    GridGenerator::hyper_cube_with_cylindrical_hole(inclusion_2, r2, sep/2, true);
    GridTools::shift(Point<2>(-sep/2, 0), inclusion_2);

    GridGenerator::merge_triangulations(inclusion_1, inclusion_2, surface);

    HyperBallBoundary<2>* inclusion_boundary_1 = new HyperBallBoundary<2>(Point<2>(sep/2, 0), r1);
    surface.set_boundary(5, *inclusion_boundary_1); 

    for(unsigned int aa = 0; aa < 10; ++aa){
      typename Triangulation<2>::active_cell_iterator cell_1 = surface.begin_active(), endc_1 = surface.end();
        for(; cell_1 != endc_1; ++cell_1)
          for(unsigned int l_1 = 0; l_1 < GeometryInfo<2>::lines_per_cell; ++l_1){
            Point<2> edge_center_1 = cell_1 -> line(l_1) -> center(true, true);
              if(cell_1 -> line(l_1) -> at_boundary()){
                if(sqrt((std::pow((edge_center_1[0]-sep/2), 2))+std::pow(edge_center_1[1], 2)) <= (r1))
                  cell_1 -> line(l_1) -> set_all_boundary_ids(5);
              }
          }
    }

    HyperBallBoundary<2>* inclusion_boundary_2 = new HyperBallBoundary<2>(Point<2>(-sep/2, 0), r2);
    surface.set_boundary(6, *inclusion_boundary_2); 
    
    for(unsigned int aa = 0; aa < 10; ++aa){
      typename Triangulation<2>::active_cell_iterator cell_2 = surface.begin_active(), endc_2 = surface.end();
        for(; cell_2 != endc_2; ++cell_2)
          for(unsigned int l_2 = 0; l_2 < GeometryInfo<2>::lines_per_cell; ++l_2){
            Point<2> edge_center_2 = cell_2 -> line(l_2) -> center(true, true);
              if(cell_2 -> line(l_2) -> at_boundary()){
                if(sqrt((std::pow((edge_center_2[0]+sep/2), 2))+std::pow(edge_center_2[1], 2)) <= (r2))
                  cell_2 -> line(l_2) -> set_all_boundary_ids(6);
              }
          }
    }

    typename Triangulation<2>::active_cell_iterator 
    cell = surface.begin_active(), 
    endc = surface.end();
    for (; cell != endc; ++cell){
      for(unsigned int i = 0; i < GeometryInfo<2>::vertices_per_cell; ++i){
        Point<2> &v = cell -> vertex(i);
        if(v(0) == sep)
          v(0) = x/2; 
        if(v(0) == -sep)
          v(0) = -x/2;
        if(v(1) == sep/2)
          v(1) = y/2; 
        if(v(1) == -sep/2)
          v(1) = -y/2;              
      
      }

    }

    if(first_run == false){ 
      surface.clear();
      delete inclusion_boundary_1;
      delete inclusion_boundary_2;
    }

}
    