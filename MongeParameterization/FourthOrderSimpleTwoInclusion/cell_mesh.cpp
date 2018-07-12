#include "fourthorder.h"

void rectangle_with_cylindrical_hole(Triangulation<2> &triangulation,
                                     double sep, 
                                     const double inner_radius,
                                     const double x,
                                     const double y,
                                     const int n,
                                     bool colorize)
   {
    triangulation.clear();
     const int dim = 2;
 
     Assert(inner_radius < x,
            ExcMessage("outer_radius has to be bigger than inner_radius."));
      Assert(inner_radius < y,
            ExcMessage("outer_radius has to be bigger than inner_radius."));
 
     Point<dim> center;
     // We create an hyper_shell in two dimensions, and then we modify it.
     GridGenerator::hyper_shell(triangulation, center, inner_radius, y, 8);
     Triangulation<dim>::active_cell_iterator
     cell = triangulation.begin_active(),
     endc = triangulation.end();
     std::vector<bool> treated_vertices(triangulation.n_vertices(), false);
     for (; cell != endc; ++cell)
       {
         for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
           if (cell->face(f)->at_boundary())
             {
               for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v)
                 {
                   unsigned int vv = cell->face(f)->vertex_index(v);
                   if (treated_vertices[vv] == false)
                     {
                      if(n == 1){
                        treated_vertices[vv] = true;
                        switch (vv)
                          {
                          case 1:
                            cell->face(f)->vertex(v) = center+Point<dim>(x/2, y/2);
                            break;
                          case 3:
                            cell->face(f)->vertex(v) = center+Point<dim>(-sep/2, y/2);
                            break;
                          case 5:
                            cell->face(f)->vertex(v) = center+Point<dim>(-sep/2, -y/2);
                            break;
                          case 7:
                            cell->face(f)->vertex(v) = center+Point<dim>(x/2, -y/2);
                          default:
                            break;

                          }
                      }  
                      else{
                        treated_vertices[vv] = true;
                        switch (vv)
                          {
                          case 1:
                            cell->face(f)->vertex(v) = center+Point<dim>(sep/2, y/2);
                            break;
                          case 3:
                            cell->face(f)->vertex(v) = center+Point<dim>(-x/2, y/2);
                            break;
                          case 5:
                            cell->face(f)->vertex(v) = center+Point<dim>(-x/2, -y/2);
                            break;
                          case 7:
                            cell->face(f)->vertex(v) = center+Point<dim>(sep/2, -y/2);
                          default:
                            break;
                          }
                      } 
                     }
                 }
             }
       }
     double epsx = 1e-3 * x;
     double epsy = 1e-3 * y;

     cell = triangulation.begin_active();
     for (; cell != endc; ++cell)
       {
         for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
           if (cell->face(f)->at_boundary())
             {
               double dx = cell->face(f)->center()(0) - center(0);
               double dy = cell->face(f)->center()(1) - center(1);
               if (colorize)
                 {
                   if (std::abs(dx + x) < epsx)
                     cell->face(f)->set_boundary_id(0);
                   else if (std::abs(dx - x) < epsx)
                     cell->face(f)->set_boundary_id(1);
                   else if (std::abs(dy + y) < epsy)
                     cell->face(f)->set_boundary_id(2);
                   else if (std::abs(dy - y) < epsy)
                     cell->face(f)->set_boundary_id(3);
                   else
                     cell->face(f)->set_boundary_id(4);
                 }
               else
                 {
                   double d = (cell->face(f)->center() - center).norm();
                   if (d-inner_radius < 0)
                     cell->face(f)->set_boundary_id(1);
                   else
                     cell->face(f)->set_boundary_id(0);
                 }
             }
       }
   }

void FourthOrder::cell_mesh(double r1, double r2, double sep, double x, double y, bool first_run){
	surface.clear();

	Triangulation<2> inclusion_1;
	rectangle_with_cylindrical_hole(inclusion_1, sep, r1, x, y, 1, false);
  GridTools::shift(Point<2>(sep/2, 0), inclusion_1);

	Triangulation<2> inclusion_2;
	rectangle_with_cylindrical_hole(inclusion_1, sep, r2, x, y, 1, false);
  GridTools::shift(Point<2>(-sep/2, 0), inclusion_2);

  GridGenerator::merge_triangulations(inclusion_1, inclusion_2, surface);
   
	std::ofstream out("twoDgrids/testgrid" + std::to_string(sep) + ".eps");
	GridOut cell_mesho;
	cell_mesho.write_eps(surface, out);

	if(first_run == false){
		surface.clear();
		//delete pointers
	}
}