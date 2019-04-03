// Class to define a vertex in 2D
//		Members: cells, edges, connected vertices...
//		Methods: index, coordinates
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include <vertex2d.hpp>
#include <mesh2d.hpp>
#include <iostream>
#include <boost/timer.hpp>

using namespace HArDCore2D;
Vertex2D::Vertex2D(size_t iV, Vector2d coords, Mesh2D *mesh)
    : _iV(iV), 
			_coords(coords), 
			_mesh(mesh),
			_cells(0),
			_edges(0),
			_vertices(0),
			_boundary(false) {
				// Do nothing
				}

Vertex2D::~Vertex2D() {}

Cell2D *Vertex2D::cell(size_t i) const {
    if (i < _cells.size()) {
        return _cells[i];
    } else {
        throw "No cell at vertex local index";
    }
}

Edge2D *Vertex2D::edge(size_t i) const {
    if (i < _edges.size()) {
        return _edges[i];
    } else {
        throw "No edge at vertex local index";
    }
}

Vertex2D *Vertex2D::vertex(size_t i) const {
    if (i < _vertices.size()) {
        return _vertices[i];
    } else {
        throw "No vertex at vertex local index";
    }
}

void Vertex2D::add_cell(Cell2D *cell) {
    _cells.push_back(cell);
}

void Vertex2D::add_edge(Edge2D *edge) {
    _edges.push_back(edge);
}

void Vertex2D::add_vertex(Vertex2D *vertex) {
    _vertices.push_back(vertex);
}

void Vertex2D::set_boundary(bool val) {
		_boundary = val;
	}

