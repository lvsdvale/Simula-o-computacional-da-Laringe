
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
#include <cmath>
#include <fstream>
#include <iostream>
#define CGAL_MESH_2_OPTIMIZER_VERBOSE
//#define CGAL_MESH_2_OPTIMIZERS_DEBUG
//#define CGAL_MESH_2_SIZING_FIELD_USE_BARYCENTRIC_COORDINATE

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

	


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
  CDT cdt;
  float H,L,f,l,h,aux=0,c;
  std::vector<float> Parametros;
  std::vector<Vertex_handle> v;
  std::string name,txtname;
  std::cout<<"digite o arquivo txt Base para a laringe"<<std::endl;
  std::cin>>txtname;
  txtname.append(".txt");
  std::fstream file;
  file.open(txtname);
  if(file.is_open()){
   while (true) {
    double value;
    if (!(file >> value)) {
      break;
    }
    Parametros.push_back(value);
   }  
  }
  else{
      std::cout<<"Nao foi possivel abrir o arquivo"<<std::endl;
      return 0;
  }

  std::cout<<"Recolhendo dados do "<<txtname<<std::endl;
  H=Parametros[0];
  std::cout<<"altura maxima: "<<H<<std::endl;
  L=Parametros[1];
  std::cout<<"tamanho horizontal maximo: "<<L<<std::endl;
  f=Parametros[2];
  std::cout<<"tamanho horinzontal até o furo: "<<f<<std::endl;
  l=Parametros[3];
  std::cout<<"tamanho horizontal do furo: "<<l<<std::endl;
if((l+f)>L){
std::cout<<"tamanhos de furo e distancia até ele invalidos altere para tamanhos validos "<<std::endl;
return 0;
}
h=Parametros[4];
std::cout<<"altura do furo:"<<h<<std::endl;
if(2*h>H){
std::cout<<"tamanho de furo invalido altere para um tamanho valido"<<std::endl;
}
c=Parametros[5];
std::cout<<"criterio de densidade de furos (0<c<1):"<<c<<std::endl;
std::cout<<"digite um nome para o arquivo VTK gerado"<<std::endl;
std::cin>>name;
name.append(".vtu");

  
 v.push_back(cdt.insert(Point(0,H)));
 v.push_back(cdt.insert(Point(0,0)));
 v.push_back(cdt.insert(Point(f,0)));
 //v.push_back(cdt.insert(Point(f,h)));
  for(int i = 180;i>90;i--){
    aux++;
    v.push_back(cdt.insert(Point(f+l+l*cos(i*(3.14159/180)),h*sin(i*(3.14159/180)))));
}
 v.push_back(cdt.insert(Point(f+l,h)));
 v.push_back(cdt.insert(Point(f+l,0)));
 v.push_back(cdt.insert(Point(L,0)));
 v.push_back(cdt.insert(Point(L,H)));
 v.push_back(cdt.insert(Point(f+l,H)));
 v.push_back(cdt.insert(Point(f+l,H-h)));
 for(int k = 270;k>180;k--){
    aux++;
    v.push_back(cdt.insert(Point(f+l+l*cos(k*(3.14159/180)),H+h*sin(k*(3.14159/180)))));
}

 //v.push_back(cdt.insert(Point(f,H-h)));
 v.push_back(cdt.insert(Point(f,H)));
 


for(int j=0;j<(aux+10);j++){
	if(j!=aux+9){
		cdt.insert_constraint(v[j],v[j+1]);
	}
	else{
            	cdt.insert_constraint(v[j],v[0]);
	}
  
}

  
  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout << "Meshing the triangulation..." << std::endl;
  CGAL::refine_Delaunay_mesh_2(cdt,Criteria(0.125,c));
  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout << "number of faces: " << cdt.number_of_faces() << std::endl;

  write_cdt_to_vtk_xml_file(cdt,name);
  //write_cdt_to_vtk_xml_file(cdt,txtname);
 
 }
















