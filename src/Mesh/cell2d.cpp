// Class to define a cell
//		Members: vertices, edges, neighbouring cells...
//		Methods: index, diameter, area, center of mass...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include "cell2d.hpp"
#include "mesh2d.hpp"
#include <cmath>
#include <math.h>
#include <iostream>
#include <utils.hpp>
#include <boost/timer.hpp>

#include <stdlib.h>     /* exit, EXIT_FAILURE */

using namespace HArDCore2D;

Cell2D::Cell2D(size_t iC, std::vector<size_t> vertex_ids, Mesh2D *mesh)
    : _iC(iC),
      _vertex_ids(vertex_ids),
      _mesh(mesh),
	 		_edges(0),
	 		_vertices(0),
	  	_neighbours(0),
      _boundary(false) {
		    // need to initialise the cells edge - while we're at it why not calculate
		    // the edge midpoints
		    std::vector<size_t> tmp = {1,1};
		    for (size_t i = 0; i < vertex_ids.size(); i++) {
					// grab vertex i and i+1 (modulo nb of vertices in cell)
	        size_t j = i + 1;
	        if (j >= vertex_ids.size()) {
	            j = 0;
	        }
					// Vertices i,j appear counter-clockwise in cell
	        tmp[0] = vertex_ids[i];
	        tmp[1] = vertex_ids[j];

	        Edge2D *edge = _mesh->add_edge(tmp, this);
	        _edges.push_back(edge);

					// Add vertex to cell
					Vertex2D *vertex = _mesh->vertex(vertex_ids[i]);
					_vertices.push_back(vertex);

    		}	
				calc_cell_geometry_factors();

			}

Cell2D::~Cell2D() {}

std::vector<Edge2D *> Cell2D::get_edges() const { return _edges; }

std::vector<Vertex2D *> Cell2D::get_vertices() const { return _vertices; }

std::vector<Cell2D *> Cell2D::get_neighbours() const { return _neighbours; }

Edge2D *Cell2D::edge(size_t i) const {
    if (i < _edges.size()) {
        return _edges[i];
    } else {
        throw "No edge at local index";
    }
}

Vertex2D *Cell2D::vertex(size_t i) const {
    if (i < _vertices.size()) {
        return _vertices[i];
    } else {
        throw "No vertex at local index";
    }
}

Cell2D *Cell2D::neighbour(size_t i) const {
    if (i < _neighbours.size()) {
        return _neighbours[i];
    } else {
        throw "No neighbour at local index";
    }
}

Vector2d Cell2D::edge_normal(size_t i) {
    //because we need to know the order of the indices in the cell we should just recreate the edge indicies here!
    size_t k = i;
    size_t j = i+1;
    if (j>=n_vertices()){
        j = 0;
        }
    Vector2d v1 = vertex(k)->coords();
    Vector2d v2 = vertex(j)->coords();

    Vector2d normal = Vector2d((v2 - v1).y(), -(v2 - v1).x());
    return normal.normalized();
}

Vector2d Cell2D::ari_coords() const {
    double x = 0.0;
    double y = 0.0;
    // calculate the arithmetic coordinates using the calculated edge midpoints.
    // Currently this is just computed when needed but it could be saved or even
    // precomputed when the midpoints are calculated
    for (size_t i = 0; i < n_edges(); i++) {
        Edge2D *e = edge(i);
        auto mp = e->center_mass();
        x += mp.x();
        y += mp.y();
    }

    x /= n_edges();
    y /= n_edges();
    Vector2d ari(x, y);
    return ari;
}

bool Cell2D::is_neighbour(const Cell2D *rhs) const {
    // iterate over the nodes indices in the cells they are shared then
    // these
    // are neighbours
    size_t n = 0;
    if (rhs->global_index() == global_index()) {
        return false;
    }
    for (size_t i = 0; i < n_vertices(); i++) {
        for (size_t j = 0; j < rhs->n_vertices(); j++) {
            if (vertex(i) == rhs->vertex(j)) {
                n++;  // return true;
            }
        }
    }
    if (n >= 2) {
        return true;
    }
    return false;
}

bool Cell2D::add_neighbour(Cell2D *neigh) {
    _neighbours.push_back(neigh);
    return true;
}

bool Cell2D::calc_cell_geometry_factors() {
    std::vector<Vertex2D *> vlist = get_vertices();
    size_t nFV = vlist.size();
    _cell_diam = 0.0;
    for (size_t iVl = 0; iVl < nFV; iVl++) {
				Vector2d p1 = vlist[iVl]->coords();
        for (size_t jVl = iVl + 1; jVl < nFV; jVl++) {
						Vector2d p2 = vlist[jVl]->coords();
            _cell_diam = std::max(_cell_diam, (p1 - p2).norm());
        }
    }
		Vector2d xc = Vector2d::Zero();
    double de(0.), ar(0.);

    if (nFV == 3) {  // cell is a triangle
				Vector2d v0 = vlist[1]->coords() - vlist[0]->coords();
				Vector2d v1 = vlist[2]->coords() - vlist[0]->coords();
				de = v0(0) * v1(1) - v0(1) * v1(0);
        ar = sqrt(de * de);
        // be careful that ar is twice the triangle area!
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
					xc += vlist[ilV]->coords() / 3.;
        }
				_center_mass = xc;

    } else {  // F is a planar polygon
        // get sub-triangle vertices
        Vector2d ari_center = ari_coords();
        for (size_t ilV = 0; ilV < nFV; ++ilV) {
						size_t ilVnext = (ilV + 1) % nFV;
						Vector2d v1 = vlist[ilV]->coords();
						Vector2d v2 = vlist[ilVnext]->coords();
						Vector2d x = v1 - ari_center;
						Vector2d y = v2 - ari_center;

            double tmp_de = x(0) * y(1) - x(1) * y(0);
            double tmp_ar = sqrt(tmp_de * tmp_de);
            ar += tmp_ar;
            // be careful that: tmp_ar is twice the triangle area!
            // (the coefficient should be sub_area_triangle/6)
						xc += (ari_center + v1 + v2) * tmp_ar / 3;
        }
				_center_mass = xc / ar;
    }
    // --------------------------
    _cell_area = ar / 2.;

    return true;  
}

void Cell2D::set_boundary(bool val) {
		_boundary = val;
	}


size_t Cell2D::shared_edge(size_t i) {
    Edge2D *e = edge(i);
		size_t neighbour_cell = _mesh->n_cells() + 1;
    if (e->is_boundary()) {
        neighbour_cell = _mesh->n_cells();
    }
    //get the index for both cells that belong to this edge
    size_t iC1 = e->cell(0)->global_index();
    size_t iC2 = e->cell(1)->global_index();
    //check which one isn't the current cell 
    if (iC1 != global_index()) {
        neighbour_cell = iC1;
    }
    if (iC2 != global_index()) {
        neighbour_cell = iC2;
    }

		return neighbour_cell;
}


