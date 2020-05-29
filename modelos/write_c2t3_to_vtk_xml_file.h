#ifndef Mesh_3_example_write_c3t3_to_vtk_xml_file_h
#define Mesh_3_example_write_c3t3_to_vtk_xml_file_h

#include<fstream>

namespace CGAL {
    template<class C2t3>
    bool write_c2t3_to_vtk_xml_file(const C2t3 &c2t3, const std::string &file_name)
    {
        typedef typename C2t3::Triangulation Tr;
        typedef typename Tr::Finite_facets_iterator Cell_iterator;
        typedef typename Tr::Finite_vertices_iterator Vertex_iterator;
        
        // Domain
        typedef Exact_predicates_inexact_constructions_kernel K;
        typedef K::FT FT;
        typedef K::Point_3 Point;
        
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
        Tr triangulation = c2t3.triangulation();
        int num_vertices = triangulation.number_of_vertices();
        int num_cells = c2t3.number_of_facets();
        
        indent += indent_unit;
        vtk_file << indent + "<Piece NumberOfPoints=\"" << num_vertices << "\" ";
        vtk_file << "NumberOfCells=\"" << num_cells << "\">" << std::endl;
        
        // Write vertices
        indent += indent_unit;
        vtk_file << indent + "<Points>" << std::endl;
        
        indent += indent_unit;
        vtk_file << indent;
        vtk_file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
        
        std::map<Point, int> V;
        int i=0;
        indent += indent_unit;
        
        for (Vertex_iterator it=triangulation.finite_vertices_begin(); it != triangulation.finite_vertices_end(); ++it)
        {
            vtk_file << indent;
            vtk_file << it->point().x() << " " << it->point().y() << " " << it->point().z() << std::endl;
            V[it->point()] = i;
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
        Cell_iterator it;
        for (it = triangulation.finite_facets_begin(); it != triangulation.finite_facets_end(); ++it)
        {
            if (it->first->is_facet_on_surface(it->second)) {
                vtk_file << indent;
                vtk_file << V[it->first->vertex((it->second+1)%4)->point()] << " ";
                vtk_file << V[it->first->vertex((it->second+2)%4)->point()] << " ";
                vtk_file << V[it->first->vertex((it->second+3)%4)->point()] << " " << std::endl;
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
    
} // end namespace CGAL

#endif
