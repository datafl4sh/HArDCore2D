// Creates a .vtu file of the solution, to be visualised by paraview
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include <mesh2d.hpp>
#include <stdio.h>
#include <utils.hpp>
#include <vtu_writer2d.hpp>
using namespace HArDCore2D;

VtuWriter2D::VtuWriter2D(Mesh2D* mesh) : _mesh(mesh) {}
void VtuWriter2D::write_header(FILE* pFile) {
    fprintf(pFile, "%s", "<VTKFile type=\"UnstructuredGrid\"");
    fprintf(pFile, "%s", " version=\"0.1\"");
    fprintf(pFile, "%s\n", " byte_order=\"LittleEndian\">");
    fprintf(pFile, "%s\n", "   <UnstructuredGrid>");
    fprintf(pFile, "%s", "      <Piece  ");
    fprintf(pFile, "%s%i%s", "NumberOfPoints=\"", (int)_mesh->n_vertices(), "\"  ");
    fprintf(pFile, "%s%i%s\n", "NumberOfCells=\"", (int)_mesh->n_cells(), "\">");
}

void VtuWriter2D::write_vertices(FILE* pFile) {
    // VERTICES
    fprintf(pFile, "%s\n", "         <Points>");
    fprintf(pFile, "%s%s%s\n",
            "            <DataArray type=\"Float32\" NumberOfComponents=\"",
            "3", "\" format=\"ascii\">");
    for (size_t iV = 0; iV < _mesh->n_vertices(); iV++) {
        // vertices coords (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1}
        // z_{i=1}
        //....
        // x_{i=n} y_{i=n} z_{i=n}) i is the node index
        auto v = _mesh->vertex(iV)->coords();
        fprintf(pFile, "%e %e %e", v(0), v(1), 0.0);
        fprintf(pFile, "%s\n", "");
    }

    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "         </Points>");
}

void VtuWriter2D::write_vertices(FILE* pFile, Eigen::VectorXd data) {
    // VERTICES
    fprintf(pFile, "%s\n", "         <Points>");
    fprintf(pFile, "%s%s%s\n",
            "            <DataArray type=\"Float32\" NumberOfComponents=\"",
            "3", "\" format=\"ascii\">");
    for (size_t iV = 0; iV < _mesh->n_vertices(); iV++) {
        // vertices coords (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1}
        // z_{i=1}
        //....
        // x_{i=n} y_{i=n} z_{i=n}) i is the node index
        auto v = _mesh->vertex(iV)->coords();
        fprintf(pFile, "%e %e %e", v(0), v(1), data(iV));
        fprintf(pFile, "%s\n", "");
    }

    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "         </Points>");
}


void VtuWriter2D::write_cells(FILE* pFile) {
    fprintf(pFile, "%s\n", "         <Cells>");
    fprintf(pFile, "%s\n",
            "            <DataArray type=\"Int32\" Name=\"connectivity\" "
            "format=\"ascii\">");
    // Lists all vertices of all nodes
    for (size_t iC = 0; iC < _mesh->n_cells(); iC++) {
        std::vector<Vertex2D *> vertices = _mesh->cell(iC)->get_vertices();
        for (size_t iV = 0; iV < vertices.size(); iV++) {
            fprintf(pFile, "%i ", (int)vertices[iV]->global_index());
        }
        fprintf(pFile, "%s\n", "");

    }
    fprintf(pFile, "%s\n", "");
    fprintf(pFile, "%s\n", "            </DataArray>");

    fprintf(pFile, "%s\n",
            "            <DataArray type=\"Int32\" Name=\"offsets\" "
            "format=\"ascii\">");
    size_t prev = 0;
    // probably could be done more efficiently like building fofsets
    // array when
    // making cells list but...
    for (size_t iT = 0; iT < _mesh->n_cells(); iT++) {
        prev += _mesh->cell(iT)->n_vertices();  // offsets are the end position of the cell in the nodes array?

        fprintf(pFile, "%i ", (int)prev);
    }
    fprintf(pFile, "%s\n", "");

    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n",
            "            <DataArray type=\"Int32\" Name=\"types\" "
            "format=\"ascii\">");
    // probably could be done more efficiently like building fofsets
    // array when
    // making cells list but...
    for (size_t iT = 0; iT < _mesh->n_cells(); iT++) {
        fprintf(pFile, "%i ", 7);
    }
    fprintf(pFile, "%s\n", "");
    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "            </Cells>");
}

void VtuWriter2D::write_point_scalar_property(FILE* pFile,
                                                  Eigen::VectorXd data,
                                                  std::string name) {
    fprintf(pFile, "%s\n", "         <PointData Scalars=\"count\">");
    fprintf(pFile, "%s%s%s\n", "            <DataArray type=\"Float64\" Name=\"",
            name.c_str(),
            "\" "
            "format=\"ascii\">");
    for (size_t i = 0; i < size_t(data.size()); i++) {
        fprintf(pFile, "%e ", data[i]);
    }
    fprintf(pFile, "%s\n", "");

    fprintf(pFile, "%s\n", "            </DataArray>");

    fprintf(pFile, "%s\n", "         </PointData>");
}

void VtuWriter2D::write_footer(FILE* pFile) {
    fprintf(pFile, "%s\n", "      </Piece>  ");
    fprintf(pFile, "%s\n", "   </UnstructuredGrid>");
    fprintf(pFile, "%s\n", "</VTKFile>");
}
bool VtuWriter2D::write_to_vtu(std::string filename, Eigen::VectorXd sol_vertex, bool dimen) {
    FILE* pFile = fopen(filename.c_str(), "w");
    write_header(pFile);
	if (dimen) { 
		write_vertices(pFile,sol_vertex);
	} else {
		write_vertices(pFile);
	}
	write_point_scalar_property(pFile, sol_vertex, "solution");
    write_cells(pFile);
    write_footer(pFile);
    return true;
}

bool VtuWriter2D::write_to_vtu(std::string filename) {
    FILE* pFile = fopen(filename.c_str(), "w");
    write_header(pFile);
    write_vertices(pFile);
    write_cells(pFile);

    write_footer(pFile);
    return true;
}

