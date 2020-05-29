#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>
#include <CGAL/draw_triangulation_2.h>
#include <stdlib.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <string>
#define CGAL_MESH_2_OPTIMIZER_VERBOSE
//#define CGAL_MESH_2_OPTIMIZERS_DEBUG
//#define CGAL_MESH_2_SIZING_FIELD_USE_BARYCENTRIC_COORDINATE

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;
bool write_cdt_to_vtk_xml_file(CDT cdt, const std::string &file_name)    	
{
        
        typedef typename CDT::Finite_faces_iterator Cell_iterator;
        typedef typename CDT::Finite_vertices_iterator Vertex_iterator;
        
        // Domain
        typedef K::FT FT;
        typedef K::Point_3 Point3;
        
       // check that file extension is "vtu"
        CGAL_assertion(file_name.substr(file_name.length()-4,4) == ".vtu");
        
        // open file
        std::ofstream vtk_file(file_name.c_str());
        
        // header
        vtk_file << "<VTKFile type=\"UnstructuredGrid\" ";
        vtk_file << "version=\"0.1\" ";
        vtk_file << "byte_order=\"BigEndian\">" << std::endl;
        
        int indent_size = 2;
        std::string indent_unit(indent_size, ' ');
        std::string indent = indent_unit;
        vtk_file << indent + "<UnstructuredGrid>" << std::endl;
        
        // write mesh
   
        int num_vertices = cdt.number_of_vertices();
        int num_cells;
  for (Cell_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
        {
             if (it->is_in_domain()) {
                num_cells++;
               
           }
              
        }

        
        indent += indent_unit;
        vtk_file << indent + "<Piece NumberOfPoints=\"" << num_vertices << "\" ";
        vtk_file << "NumberOfCells=\"" << num_cells << "\">" << std::endl;
        
        // Write vertices
        indent += indent_unit;
        vtk_file << indent + "<Points>" << std::endl;
        
        indent += indent_unit;
        vtk_file << indent;
        vtk_file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
        
        std::map<Vertex_handle, int> V;
        int i=0;
        indent += indent_unit;
        
        for (Vertex_iterator it=cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it)
        {
            vtk_file << indent;
            vtk_file << it->point().x() << " " << it->point().y() << " " << 0 << std::endl;
            V[it.base()] = i;
            ++i;
            //std::cout<<"interação de vertices numero:"<<i<<std::endl;
        }
        
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</DataArray>" << std::endl;

        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</Points>" << std::endl;

        // Write tetrahedra
        vtk_file << indent << "<Cells>" << std::endl;
        
        indent += indent_unit;
        vtk_file << indent;
        vtk_file << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">";
        vtk_file << std::endl;
        
        indent += indent_unit;
        
        for (Cell_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
        {
             if (it->is_in_domain()) {
                vtk_file << indent;
                vtk_file << V[it->vertex(0)] << " ";
                vtk_file << V[it->vertex(1)] << " ";
                vtk_file << V[it->vertex(2)] << " " << std::endl;
               
               
           }
               
                
        }
       
        
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</DataArray>" << std::endl;
        
        // offsets
        // every element is a three node triangle so all offsets are multiples of 3
        vtk_file << indent;
        vtk_file << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">";
        vtk_file << std::endl;
        i = 3;
        indent += indent_unit;
        for (int j = 0; j < num_cells; ++j)
        {
            vtk_file << indent << i << std::endl;
            i += 3;
        }
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</DataArray>" << std::endl;
        
        // cell types (type 5 is a three node triangle)
        vtk_file << indent;
        vtk_file << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">";
        vtk_file << std::endl;
        indent += indent_unit;
        for (int j = 0; j < num_cells; ++j)
        {
            vtk_file << indent << "5" << std::endl;
        }
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</DataArray>" << std::endl;
        
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</Cells>" << std::endl;
        
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</Piece>" << std::endl;
        
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << indent + "</UnstructuredGrid>" << std::endl;
        
        indent.erase(indent.length()-indent_size, indent_size);
        vtk_file << "</VTKFile>" << std::endl;
        
        return true;
    }

int main()
{
  std::string name;
  std::cin>>name;
  CDT cdt;
  Vertex_handle va = cdt.insert(Point(2,0));
  Vertex_handle vb = cdt.insert(Point(0,2));
  Vertex_handle vc = cdt.insert(Point(-2,0));
  Vertex_handle vd = cdt.insert(Point(0,-2));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  va = cdt.insert(Point(3,3));
  vb = cdt.insert(Point(-3,3));
  vc = cdt.insert(Point(-3,-3));
  vd = cdt.insert(Point(3,0-3));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  std::list<Point> list_of_seeds;

  list_of_seeds.push_back(Point(0, 0));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing the domain..." << std::endl;
  CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                               Criteria(.125,.1),true);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
  int mesh_faces_counter = 0;
  for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
      fit != cdt.finite_faces_end(); ++fit) 
  {
    if(fit->is_in_domain()) ++mesh_faces_counter;
  }



fenics 
sha256:08e2bf44f19aeb0db5f5319f2cdb15808e4453d788ffd1971c6b6956953c0228
sha256:08e2bf44f19aeb0db5f5319f2cdb15808e4453d788ffd1971c6b6956953c0228

  std::cout << "Number of faces in the mesh domain: " << mesh_faces_counter << std::endl;
  write_cdt_to_vtk_xml_file(cdt,name);
}
