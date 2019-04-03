// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "mesh2d_builder.hpp"
#include <iostream>
#include <utils.hpp>
#include <deque>
#include <boost/timer.hpp>

using namespace HArDCore2D;
Mesh2DBuilder::Mesh2DBuilder() {}
Mesh2D* Mesh2DBuilder::build_the_mesh(
    std::vector<std::vector<double> > vertices,
    std::vector<std::vector<size_t> > cells) {
    if (vertices.size() > 0 && cells.size() > 0) {
        _mesh = new Mesh2D;  // make a pointer to the mesh so that it outlives
                             // the builder
        std::cout << "Mesh: ";

				// Create vertices
        size_t iG = 0;
        for (std::vector<std::vector<double> >::iterator it = vertices.begin(); it != vertices.end(); it++) {
            Vertex2D* vertex = new Vertex2D(iG, Vector2d((*it)[0], (*it)[1]), _mesh);
            _mesh->add_vertex(vertex);
            iG++;
        }

				// Create cells
        double total_area = 0.0;
        iG = 0;
        for (std::vector<std::vector<size_t> >::iterator it = cells.begin(); it != cells.end(); it++) {
            // first entry is the number of nodes so skip it
            std::vector<size_t> vertex_ids;
            for (size_t i = 1; i < (*it).size(); i++) {
                vertex_ids.push_back((*it)[i]);
            }
            Cell2D* cell = new Cell2D(iG, vertex_ids, _mesh);
            total_area += cell->measure();
            _mesh->add_cell(cell);
            iG++;

						// add the cell to all its vertices
						for (size_t ilV = 0; ilV < vertex_ids.size(); ilV++){
							size_t iV = vertex_ids[ilV];
							_mesh->vertex(iV)->add_cell(cell);
						}
        }

				// build boundary
				build_boundary();

        std::cout << "added " << _mesh->n_cells() << " cells; Total area= " << total_area << std::endl;
        return _mesh;
    } else {
        if (vertices.size() <= 0) {
            throw "Can't build mesh vertices is empty. Check the input file";
        }
        if (cells.size() <= 0) {
            throw "Can't build mesh cells is empty. Check the input file";
        }
    }
    return NULL;
}

void Mesh2DBuilder::build_boundary() {
		// When the mesh is built, boundary edges are already identified (see Edge2D::add_cell)
		// Here we fill in the _boundary variables of the cells and vertices, and the lists of boundary
		// edges, cells and vertices
		for (size_t iC = 0; iC < _mesh->n_cells(); iC++) {
			Cell2D* cell = _mesh->cell(iC);
			for (size_t ilE = 0; ilE < cell->n_edges(); ilE++) {
				Edge2D* edge = cell->edge(ilE);
				if (edge->is_boundary()) {
					// The cell has a boundary edge, so it is a boundary cell
					cell->set_boundary(true);
					_mesh->add_b_cell(cell);
					// We also add the edge to the boundary edges
					_mesh->add_b_edge(edge);
				}
			}
			// If we have a boundary cell, we explore its vertices and those connected to 
			// boundary edges are boundary vertices
			if (cell->is_boundary()) {
				for (size_t ilV = 0; ilV < cell->n_vertices(); ilV++)	{
					Vertex2D* vertex = cell->vertex(ilV);
					for (size_t i = 0; i < 2; i++){
						if (vertex->edge(i)->is_boundary()){
							vertex->set_boundary(true);
							_mesh->add_b_vertex(vertex);
						}
					}
				}

			}
		}
	
}


