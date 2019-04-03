// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef MESH2D_BUILDER_HPP
#define MESH2D_BUILDER_HPP
#include "mesh2d.hpp"
#include "cell2d.hpp"
#include "edge2d.hpp"
#include "vertex2d.hpp"
#include <string>
#include <vector>
#include <Eigen/Dense>

#include <stdlib.h>     /* exit, EXIT_FAILURE */

namespace HArDCore2D {

/*!
*	@addtogroup Mesh2D
* @{
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The Mesh2DBuilder class provides build tools to create a full mesh with all connectivities
class Mesh2DBuilder {
public:
    /**
    * Constructor for Mesh2DBuilder.
    */
    Mesh2DBuilder();
    /**
    * Build a Mesh2D from vertices and cells
    *
    * @param vertices vector containing the coordinates of the vertices ordered
    *using the global ordering. Note that indexes start at 0 in c++
    * @param cells vector containing vectors with  the global indexes of
    *vertices making up a cell
    *
    * @return a pointer to the mesh that was build
    */
    Mesh2D* build_the_mesh(std::vector<std::vector<double> > vertices,
                           std::vector<std::vector<size_t> > cells);  ///< construct the connectivity in the mesh
    inline Mesh2D* mesh();  ///< getter for the mesh that was built by this builder

private:
    void build_boundary();  ///< identifies boundary cells and vertices, and compute lists of boundary cells, edges and vertices

    Mesh2D* _mesh;  ///< the mesh that is built by this builder
};
Mesh2D* Mesh2DBuilder::mesh() { return _mesh; }

/*@}*/
}
#endif /* MESH2D_BUILDER_HPP */

