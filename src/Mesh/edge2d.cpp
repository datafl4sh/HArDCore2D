// Class to define an edge
//		Members: cells, vertices...
//		Methods: index, length, midpoint...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include <edge2d.hpp>
#include <mesh2d.hpp>
#include <iostream>
#include <boost/timer.hpp>

using namespace HArDCore2D;
Edge2D::Edge2D(size_t iE, std::vector<size_t> vertices, Mesh2D *mesh, Cell2D *cell)
    : _iE(iE), 
			_vertex_ids(vertices), 
			_mesh(mesh), 
			_boundary(true),
			_cells(0) {
		    _cells.push_back(cell);
				_vertices.push_back(_mesh->vertex(vertices[0]));
				_vertices.push_back(_mesh->vertex(vertices[1]));
		    Vector2d p1 = _mesh->vertex(vertices[0])->coords();
		    Vector2d p2 = _mesh->vertex(vertices[1])->coords();

		    if (_cells[0] == NULL) {
	        std::cout << "cell nil" << std::endl;
		    }
		    _line = p1 - p2;
				_mp = (p1 + p2)/2;
}

Edge2D::~Edge2D() {}

Cell2D *Edge2D::cell(size_t i) const {
    if (i < _cells.size()) {
        return _cells[i];
    } else {
        throw "No cell at edge local index";
    }
}

Vertex2D *Edge2D::vertex(size_t i) const {
    if (i < _vertices.size()) {
        return _vertices[i];
    } else {
        throw "No vertex at edge local index";
    }
}

std::vector<Cell2D *> Edge2D::get_cells() const { return _cells; }

std::vector<Vertex2D *> Edge2D::get_vertices() const { return _vertices; }

double Edge2D::measure() const { return _line.norm(); }

Vector2d Edge2D::center_mass() const { return _mp; }

Vector2d Edge2D::tangent() const { return _line; }

void Edge2D::add_cell(Cell2D *cell) {
    if (_cells.size() < 2) {
        _cells.push_back(cell);
        // if the edge has 2 cells its not a boundary
        _boundary = false;
    } else {
        throw "An edge cannot have three cells!";
    }
}

