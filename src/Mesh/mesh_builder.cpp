// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "mesh_builder.hpp"
#include <iostream>
#include <utils.hpp>
#include <deque>

using namespace HArDCore2D;
MeshBuilder::MeshBuilder() {}
Mesh* MeshBuilder::build_the_mesh(
    std::vector<std::vector<double> > vertices,
    std::vector<std::vector<size_t> > cells) {
    if (vertices.size() > 0 && cells.size() > 0) {
        _mesh = new Mesh;  // make a pointer to the mesh so that it outlives the builder
        std::cout << "Mesh: ";

				// Create vertices
        size_t iG = 0;
        for (auto& v : vertices) {
            Vertex* vertex = new Vertex(iG, Vector2d(v[0], v[1]), _mesh);
            _mesh->add_vertex(vertex);
            iG++;
        }

				// Create cells
        double total_area = 0.0;
        iG = 0;
        for (auto& c : cells) {
            // first entry is the number of nodes so skip it
            std::vector<size_t> vertex_ids;
            for (size_t i = 1; i < c.size(); i++) {
                vertex_ids.push_back(c[i]);
            }
            Cell* cell = new Cell(iG, vertex_ids, _mesh);
            total_area += cell->measure();
            _mesh->add_cell(cell);
            iG++;

						// add the cell to all its vertices
						for (auto& vid : vertex_ids){
							_mesh->vertex(vid)->add_cell(cell);
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

void MeshBuilder::build_boundary() {
		// When the mesh is built, boundary edges are already identified (see Edge::add_cell)
		// Here we fill in the _boundary variables of the cells and vertices, and the lists of boundary
		// edges, cells and vertices
		for (auto& cell : _mesh->get_cells()) {
			for (auto& edge : cell->get_edges()) {
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
				for (auto& vertex : cell->get_vertices())	{
					for (size_t i = 0; i < 2; i++){
						if (vertex->edge(i)->is_boundary()){
							vertex->set_boundary(true);
							_mesh->add_b_vertex(vertex);
						}
					}
				}

			}
		}

		// Pass to fill in interior elements
		for (auto& cell : _mesh->get_cells()){
			if ( !(cell->is_boundary()) ){
				_mesh->add_i_cell(cell);
			}	
		}
		for (auto& edge : _mesh->get_edges()){
			if ( !(edge->is_boundary()) ){
				_mesh->add_i_edge(edge);
			}	
		}
		for (auto& vertex : _mesh->get_vertices()){
			if ( !(vertex->is_boundary()) ){
				_mesh->add_i_vertex(vertex);
			}	
		}
}


