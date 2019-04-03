// Class to define a vertex in 2D
//		Members: cells, edges, connected vertices...
//		Methods: index, coordinates
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#ifndef VERTEX2D_HPP
#define VERTEX2D_HPP
#include <cell2d.hpp>
#include <edge2d.hpp>
#include <mesh2d.hpp>
namespace HArDCore2D {  // forward declaration
class Mesh2D;
}


namespace HArDCore2D {

/*!
*	@addtogroup Mesh2D
* @{
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The Vertex2D class provides description of a vertex
class Vertex2D {
    /**
    * A class representing a vertex of a cell for a 2D mesh. Contains a
    * pointer to the mesh, the vertex coordinates, and cells, edges and vertices linked to that vertex.
    */
public:
    /**
    * Default constructor
    *
    * @param iV global vertex number
    * @param coord  coordinates of the vertex
    * @param mesh pointer to the mesh
    * @param cell pointer to the cells the vertex belongs to
    * @param edge pointer to the edges the vertex belongs to
    * @param vertex pointer to the vertex linked to this vertex by an edge
    */
    Vertex2D(size_t iV, Vector2d coords, Mesh2D *mesh);
    ~Vertex2D(); // destructor, nothing special

    inline size_t global_index() const;  ///< returns the edges global index
    inline Vector2d coords() const;			///< returns the coordinates of the vertex

    inline size_t n_cells() const;			///< returns the number of cells that contain the vertex
    inline size_t n_edges() const;			///< returns the number of edges that contain the vertex
    inline size_t n_vertices() const;			///< returns the number of vertices connected by an edge to the vertex
    inline std::vector<Cell2D *> get_cells() const;			///< returns the list of cells containing the vertex			
    inline std::vector<Edge2D *> get_edges() const;			///< returns the list of edges containint the vertex			
    inline std::vector<Vertex2D *> get_vertices() const;	///< returns the list of vertices linked to the vertex
		Cell2D *cell(size_t i) const; 			///< returns i-th cell containing the vertex
		Edge2D *edge(size_t i) const; 			///< returns i-th edge containing the vertex
		Vertex2D *vertex(size_t i) const; 			///< returns i-th vertex linked to the vertex

    void add_cell(Cell2D *cell);      ///< Add a new cell to the list of cells containing the vertex
    void add_edge(Edge2D *edge);      ///< Add a new edge to the list of edges containing the vertex
    void add_vertex(Vertex2D *vertex);      ///< Add a new vertex to the list of vertices connected by an edge to the vertex

		void set_boundary(bool val); ///< Set the _boundary value of the vertex to val

private:
    size_t _iV;
		Vector2d _coords;
    Mesh2D *_mesh;
    std::vector<Cell2D *> _cells;
    std::vector<Edge2D *> _edges;
    std::vector<Vertex2D *> _vertices;
    bool _boundary;
	
};

// ----------------------------------------------------------------------------
//                            Implementations
// ----------------------------------------------------------------------------

size_t Vertex2D::global_index() const { return _iV; }
Vector2d Vertex2D::coords() const { return _coords;}

size_t Vertex2D::n_cells() const { return _cells.size();}
size_t Vertex2D::n_edges() const { return _edges.size();}
size_t Vertex2D::n_vertices() const { return _vertices.size();}
std::vector<Cell2D *> Vertex2D::get_cells() const { return _cells;}
std::vector<Edge2D *> Vertex2D::get_edges() const { return _edges;}
std::vector<Vertex2D *> Vertex2D::get_vertices() const { return _vertices;}

/*@}*/
}
#endif /* VERTEX2D_HPP */


