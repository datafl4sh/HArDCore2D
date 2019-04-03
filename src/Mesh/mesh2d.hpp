// Class to describe a mesh.
//		Members: cells, vertices, edges...
//		Methods: h_max, add cells and edges...
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include "cell2d.hpp"
#include <edge2d.hpp>
#include "vertex2d.hpp"
#include <Eigen/Dense>
namespace HCore2D {  // forward declaration
 class Vertex2D;
}

/*!	
*	@defgroup Mesh2D 
* @brief Classes to construct and describe a 2D mesh
*/

namespace HArDCore2D {
/**
* Class which represents a 2D mesh. Contains cells, vertices, edges.
*/

using Eigen::Vector2d;

/*!
*	@addtogroup Mesh2D
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The Mesh2D class provides description of a mesh

class Mesh2D {
public:
    /**
    * default constructor for an empty mesh
    */
    Mesh2D();

    inline void set_name(std::string name);  ///< set the name of the mesh
    inline std::string get_name();  ///< getter for the edge name
    inline size_t n_cells() const;     ///< number of cells in the mesh
    inline size_t n_edges() const;     ///< number of edges in the mesh
    inline size_t n_vertices() const;  ///< number of vertices in the mesh
		inline double h_max() const;			 ///< max of diameter of cells
    inline size_t dim() const;         ///< dimension of the mesh (2)
    size_t n_b_cells() const;           ///< number of boundary cells
    size_t n_b_edges() const;           ///< number of boundary edges
    size_t n_b_vertices() const;        ///< number of boundary vertices

    inline std::vector<Cell2D*> get_cells() const;  ///< lists the cells in the mesh.
    inline std::vector<Edge2D*> get_edges() const;  ///< lists the edges in the mesh.
    inline std::vector<Vertex2D*> get_vertices() const;  ///< lists the vertices in the mesh.
    Cell2D* cell(size_t iC) const;  ///< get a constant pointer to a cell using its global index
    Edge2D* edge(size_t iE) const;   ///< get a constant pointer to an edge using its global index
    Vertex2D* vertex(size_t iV) const;   ///< get a constant pointer to a vertex using its global index

    inline std::vector<Cell2D*> get_b_cells() const;  ///< lists the boundary cells in the mesh.
    inline std::vector<Edge2D*> get_b_edges() const;  ///< lists the boundary edges in the mesh.
    inline std::vector<Vertex2D*> get_b_vertices() const;  ///< lists the boundary vertices in the mesh.
    Cell2D* b_cell(size_t iC) const;  ///< get a constant pointer to the iC-th boundary cell
    Edge2D* b_edge(size_t iE) const;   ///< get a constant pointer to the iE-th boundary edge
    Vertex2D* b_vertex(size_t iV) const;   ///< get a constant pointer to the iV-th boundary vertex

    inline bool add_cell(Cell2D* cell);  ///<  adds a cell to the mesh
    inline bool add_vertex(Vertex2D* vertex);  ///<  adds a vertex to the mesh
    Edge2D* add_edge(std::vector<size_t> vertex_ids, Cell2D* cell);  ///< add an edge to the mesh

    inline bool add_b_cell(Cell2D* cell);  ///<  adds a boundary cell to the mesh
    inline bool add_b_edge(Edge2D* edge);  ///<  adds a boundary edge to the mesh
    inline bool add_b_vertex(Vertex2D* vertex);  ///<  adds a boundary vertex to the mesh

    inline size_t next_edge_idx();  ///< gets the next global edge index
		
private:
    std::string _mesh_name;
    size_t _next_edge_idx = 0;

    // primary data: list of cells, edges, vertices...
    std::vector<Cell2D*> _cells;
    std::vector<Edge2D*> _edges;
    std::vector<Vertex2D*> _vertices;
    std::vector<Cell2D*> _b_cells;
    std::vector<Edge2D*> _b_edges;
    std::vector<Vertex2D*> _b_vertices;
		double _h_max;
    const size_t _dim = 2;
	
};



// ----------------------------------------------------------------------------
//                            Implementations
// ----------------------------------------------------------------------------

size_t Mesh2D::n_cells() const { return _cells.size(); }
size_t Mesh2D::n_edges() const { return _edges.size(); }
size_t Mesh2D::n_vertices() const { return _vertices.size(); }
double Mesh2D::h_max() const { return _h_max; }
size_t Mesh2D::dim() const { return _dim; }
size_t Mesh2D::next_edge_idx() {
    _next_edge_idx++;
    return _next_edge_idx - 1;
}
void Mesh2D::set_name(std::string name) {
    _mesh_name = name;
    return;
}
std::string Mesh2D::get_name() { return _mesh_name; }

bool Mesh2D::add_cell(Cell2D* cell) {
		// Add the cell to the mesh and update mesh size
    _cells.push_back(cell);
		_h_max = std::max(_h_max, cell->diam());

    return true;
}

bool Mesh2D::add_b_cell(Cell2D* cell) {
    _b_cells.push_back(cell);

    return true;
}

bool Mesh2D::add_b_edge(Edge2D* edge) {
    _b_edges.push_back(edge);

    return true;
}

bool Mesh2D::add_vertex(Vertex2D* vertex) {
		// Add the vertex to the mesh
    _vertices.push_back(vertex);

    return true;
}

bool Mesh2D::add_b_vertex(Vertex2D* vertex) {
    _b_vertices.push_back(vertex);

    return true;
}

std::vector<Cell2D*> Mesh2D::get_cells() const { return _cells; }
std::vector<Edge2D*> Mesh2D::get_edges() const { return _edges; }
std::vector<Vertex2D*> Mesh2D::get_vertices() const { return _vertices; }
std::vector<Cell2D*> Mesh2D::get_b_cells() const { return _b_cells; }
std::vector<Edge2D*> Mesh2D::get_b_edges() const { return _b_edges; }
std::vector<Vertex2D*> Mesh2D::get_b_vertices() const { return _b_vertices; }

/*@}*/
}

#endif /* MESH2D_HPP */

