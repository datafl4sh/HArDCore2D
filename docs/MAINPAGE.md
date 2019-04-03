
HArD::Core provides a suite of tools, in C++, to implement numerical schemes whose unknowns are polynomials in the cells and on the edges (in 2D) or faces (in 3D). The focus is on dealing on generic polytopal meshes. This documentation addresses the 2D version of HArD::Core, but similar principles are valid for the 3D version.

\tableofcontents

* [Build instructions](#build) -- How to build the libraries and the schemes
* [The mesh structure](#mesh) -- The principles of the data structure representing the mesh, and how to load a mesh
* [The HybridCore structure](#hybridcore) -- How to use the [HybridCore2d](@ref HArDCore2D::HybridCore) structure to compute integrals on the cells and faces of the mesh, and how to use the polynomial basis functions
<!--* [Example - The Hybrid Mimetic Mixed (HMM) scheme](#hmm) -- An example implementation of the HMM scheme that illustrates the use of the mesh data and the integration routines-->

<a name="build">
\section build Build instructions
</a>

\subsection buildlib Building the libraries and the schemes

To build the libraries and implemented schemes, the minimal requirements are:

* CMake version 2.6 or above (https://cmake.org/)
* A C++ compiler that supports the C++14 standard, eg. GCC (https://gcc.gnu.org/) or Clang (https://clang.llvm.org/)
* Eigen C++ library, version 3.3 or above (http://eigen.tuxfamily.org/)
* The follwing Boost C++ libraries (http://www.boost.org/)
    * Boost filesystem
    * Boost program options

Make sure that you have the development version of boost installed. On Linux, install `libboost-dev`, `libboost-filesystem-dev` and `libboost-program-options-dev` from your package manager.

The linear systems resulting from the assembled scheme are solved using the BiCGStab implementation of Eigen. An alternative (currently commented out in the schemes' implementations) is to use the MA41 solver of the HSL library. To use this alternative you will need:

* BLAS (http://www.netlib.org/blas/) and LAPACK (http://www.netlib.org/lapack/)
* GFortran (https://gcc.gnu.org/wiki/GFortran)
* HSL MA41 linear solver (http://www.hsl.rl.ac.uk/catalogue/ma41.html)

Once you have installed all of the required dependencies, set up the build directory and generate the build files by running the following from the repository root:

```
mkdir build
cd build
cmake ..
make
```

After this, `build/Schemes` will contain the executables (e.g. `hho-diffusion`) to run the schemes. These executable need to access the typ2 meshes, which they should naturally find if you put the `typ2_meshes` directory at the root of the project's files.


\subsection doco Building the Documentation

The mesh documentation is built with Doxygen (http://www.stack.nl/~dimitri/doxygen/). If you are reading this then somebody has already built it for you. If you modify the code and wish to rebuild the documentation, simply run `doxygen` from the root directory. The HTML version of the documentation is generated inside `documentation/html` and the LaTeX version is generated inside `documentation/latex` and can be compiled using the generated Makefile.


<a name="mesh">
\section mesh The mesh
</a>

\subsection meshpple Principles

After it is loaded, the mesh is represented by classes describing a vertex, an edge and a mesh: [Vertex2D](@ref HArDCore2D::Vertex2D), [Edge2D](@ref HArDCore2D::Edge2D), and [Cell2D](@ref HArDCore2D::Cell2D). Each of these classes contains methods to access useful information for the corresponding element, including other geometrical quantities it is related to. The mesh itself is represented by an element of the [Mesh2D](@ref HArDCore2D::Mesh2D) class with methods to access all the vertices, edges and cells (or a particular vertex, edge or cell). In this class, each cell has a unique identifier, numbered from 0.

For example, if `mesh_ptr` is a pointer to a Mesh2D class,
\code{.cpp}
Vertex2D* vertex = mesh_ptr->vertex(5);

Eigen::Vector2d vert_coord = vertex->coords()
\endcode
would store the coordinates of the fifth vertex into the Eigen vector vert_coord. As a generic rule, all geometrical vectors and matrices are Eigen classes. Vectors describing lists (e.g. of cells, of functions...) are instances of `std::vector<...>`.


Here is an example that loops over all cells, grabs all the edges of the cell, and loops over these edges to output their length. Here, `mesh_ptr` is a pointer to the mesh.

\code{.cpp}
// Loop over all cells of the mesh
for (size_t iC = 0; iC < mesh_ptr->n_cells() iC++) {

	// We grab the edges of the iC-th cell
	std::vector<Edge2D *> edges = mesh_ptr->cell(iC)->get_edges();

	// Loop over the edges of the cell
	for (size_t ilE = 0; ilE < cell->n_edges(); ilE++) {

		// Write the edge length on the standard output
		std::cout << "The length of edge " << ilE+1 << " in cell " << iC+1 << " is: " << edges(ilE)->measure() << "\n";
	}

}
\endcode

The mesh classes and other auxilliary classes are located inside the namespace [HArDCore2D](@ref HArDCore2D).

There is no direct access from a high-level geometrical entity to elements purely associated with lower-level entities. For example, if `mesh_ptr` is a pointer to the mesh, there is no direct method to access the coordinates of the i-th vertex of the mesh (no `mesh_ptr->coords_vertex()` exists). Instead, this is done through `mesh_ptr->vertex(i)->coords()`. This choice is deliberate as it preserves the logical organisation of the data structure, and facilitates the memorisation of the various methods. Of course, writing a wrapper providing a direct access is easy...


\subsection loading_mesh Loading a mesh

HArDCore2D can currently read meshes in the `typ2` format designed for the FVCA5 Benchmark. A short documentation describing this format is provided in the `typ2_meshes` directory (see README.pdf there). Several meshes can also be found in this directory.

A mesh file must be read using an instance of the `Mesh2DReaderTyp2` class, and then built using `Mesh2DBuilder`.  A working example is given below (assuming the executable will be in `build/Schemes` for example).

\code{.cpp}
#include "mesh2d.hpp"
#include "import_mesh2d.hpp"
#include "mesh2d_builder.hpp"

using namespace HArDCore2D;

int main() {

	// Mesh file to read
	std::string mesh_file = "../../typ2_meshes/cart5x5.typ2";

	// Read the mesh file
	Mesh2DReaderTyp2 mesh(mesh_file);
	std::vector<std::vector<double> > vertices;
	std::vector<std::vector<size_t> > cells;
	std::vector<std::vector<double> > centers;
	if (mesh.read_mesh(vertices, cells, centers) == false) {
		std::cout << "Could not open file" << std::endl;
		return false;
	};
	// Build the mesh
	Mesh2DBuilder* builder = new Mesh2DBuilder();
	Mesh2D* mesh_ptr = builder->build_the_mesh(vertices, cells);

	std::cout << "There are " << mesh_ptr->n_cells() << " cells in the mesh.\n";
}
\endcode


<i>Note</i>: the `typ2` format allows for meshes with very generic polygonal cells, including non-convex cells.
However, the builder assumes that each cell is star-shaped with respect to the isobarycenter of its vertices -- otherwise,
the calculation of the center of mass may be incorrect. Similarly, the quadrature rules (see [Quadrature rules](#quad_rules)) assume that each cell is star-shaped with respect to its center of mass.


<a name="hybridcore">
\section hybridcore The HybridCore structure
</a>

The [HybridCore](@ref HArDCore2D::HybridCore) structure encapsulates routines to create bases of polynomial spaces in each cell and on each edge, to integrate functions on these mesh entities, and to evaluate functions defined through their coefficients on the cell and edge basis functions. Start by including the structure as follows:

\code{.cpp}
#include "hybridcore2d.hpp"
\endcode

\subsection init_core Initialising the HybridCore structure

To initialise a [HybridCore](@ref HArDCore2D::HybridCore) structure, we must first have a mesh loaded as a pointer (see [loading a mesh](#loading_mesh)). The structure is then created by specifying the desired degree of the polynomial spaces on the edges and in the cells.

\code{.cpp}
// Create the HybridCore structure with polynomial degree K on the edges, and L in the cells
  HybridCore hho(mesh_ptr, K, L);
\endcode

<a name="basis">
\subsection basis Polynomial basis functions
</a>

Initialising the [HybridCore](@ref HArDCore2D::HybridCore) structure constructs basis functions for the polynomial spaces on the edges (up to degree \f$K\f$) and in the cells (up to degree \f$(K+1)\f$). The choice of the degree in the cells corresponds to the needs of certain high-order methods, such as the HHO method that requires the reconstruction of a polynomial of degree \f$(K+1)\f$ in each cell. It is also assumed that \f$L\le K+1\f$.

The basis functions are accessed through the methods [cell_basis(iT,i)](@ref HArDCore2D::HybridCore#cell_basis) and [edge_basis(iE,i)](@ref HArDCore2D::HybridCore#edge_basis) which return the i-th basis function on the iT-th cell or iE-th edge. The cell gradients are available from [cell_gradient(iT,i)](@ref HArDCore2D::HybridCore#cell_gradient); these gradients are indexed to correspond with their basis functions, which means that the first gradient will always identically be the zero vector, since it corresponds to the constant basis function.

\code{.cpp}
    const auto &phi_i = cell_basis(iT, i);
    const auto &dphi_i = cell_gradient(iT, i);	// dphi_i is the gradient of phi_i
\endcode

The basis functions are hierarchical, which means that they are constructed by increasing degree. A basis of the space of polynomials of degree \f$K\f$ in the cell is thus obtained by selecting the first \f$(K+1)(K+2)/2\f$ cell basis functions.

When a scheme has polynomial unknowns of degree \f$K\f$ on the edges and \f$L\f$ in the cells, these unknowns can be represented as vectors of coefficients on the basis functions, for example by listing all the coefficients on the basis functions in the first cell, then all the coefficients on the basis functions in the second cell, etc., and then listing all the coefficients on the basis functions in the first edge, etc. This is the choice adopted in [HHO-diffusion](@ref HArDCore2D::HHO-diffusion).

A number of convenient quantities relating to the basis functions are available in the [HybridCore](@ref HArDCore2D::HybridCore) structure follows.

Symbol name  | Meaning
------------------- | -
[nlocal_cell_dofs](@ref HArDCore2D::HybridCore#nlocal_cell_dofs)  | The number of degrees of freedom of a cell polynomial (ie. the dimension of the space of polynomials of degree \f$\le L\f$ in two variables)
[nlocal_face_dofs](@ref HArDCore2D::HybridCore#nlocal_edge_dofs)  | The number of degrees of freedom of a face polynomial (ie. the dimension of the space of polynomials of degree \f$\le K\f$ in one variable)
[nhighorder_dofs](@ref HArDCore2D::HybridCore#nhighorder_dofs) | The number of degrees of freedom of a degree \f$(K+1)\f$ cell polynomial (ie. the dimension of the space of polynomials of degree \f$\le K+1\f$ in two variables)
[ngradient_dofs](@ref HArDCore2D::HybridCore#ngradient_dofs) | The dimension of the gradient space of degree \f$(K+1)\f$ cell polynomials
[ntotal_cell_dofs](@ref HArDCore2D::HybridCore#ntotal_cell_dofs) | The total number of degrees of freedom over all cell polynomials over the entire mesh (i.e. the number of cells times the dimension of the space of polynomials of degree \f$\le L\f$ in two variables)
[ntotal_face_dofs](@ref HArDCore2D::HybridCore#ntotal_edge_dofs) | The total number of degrees of freedom over all face polynomials over the entire mesh (i.e. the number of edges times the dimension of the space of polynomials of degree \f$\le K\f$ in one variables)
[ninternal_face_dofs](@ref HArDCore2D::HybridCore#ninternal_edge_dofs) | The total number of degrees of freedom over all internal faces over the entire mesh
[ntotal_dofs](@ref HArDCore2D::HybridCore#ntotal_dofs) | The total number of cell and face degrees of freedom over the entire mesh


\subsection integrate_mesh Integration over cells and edges

The [HybridCore](@ref HArDCore2D::HybridCore) structure provides routines to integrate generic functions on the cells and the edges. These routines are however expensive as they re-compute the quadrature points and weights every time they are called. They should therefore only be used with parsimony; computing quadrature nodes and values of basis functions at these nodes is more efficient, see [Quadratures rules](#quad_rules).

For example, to integrate \f$f(x,y) = x^2 + y^2 \f$ over the cell number iT:

\code{.cpp}
auto integral = integrate_over_cell(iT, [](auto x, auto y) {
        return x*x+y*y;
      });
\endcode

Basis functions can also be integrated:

\code{.cpp}
const auto& phi_i = edge_basis(iT, i);
const auto& phi_j = cell_basis(iT, j);
auto integral = hho.integrate_over_face(iT, [&phi_i, &phi_j](auto x, auto y, auto z) {
  return phi_i(x,y) * phi_j(x,y);
});
\endcode


<!--To illustrate, the following code computes a mass matrix of the cell basis for a particular region `iR`. We use the [Eigen](http://eigen.tuxfamily.org/) library for the matrix.-->

<!--\code{.cpp}-->
<!--Eigen::MatrixXd M(hho.nlocal_cell_dofs, hho.nlocal_cell_dofs);-->

<!--for (size_t i = 0; i < hho.nlocal_cell_dofs; i++) {-->
<!--  const auto& phi_i = hho.cell_basis(iR, i);-->
<!--  for (size_t j = 0; j < hho.nlocal_cell_dofs; j++) {-->
<!--    const auto& phi_j = hho.cell_basis(iR, j);-->
<!--    M(i, j) = hho.integrate_over_region(iR, [&phi_i, &phi_j](auto x, auto y, auto z) {-->
<!--      return phi_i(x,y,z) * phi_j(x,y,z);-->
<!--    });-->
<!--  }-->
<!--}-->
<!--\endcode-->


\subsection quad_rules Quadrature rules

HArD::Core deals with quite arbitrary cell geometries. As a consequence, no reference element can be used, and the quadrature rules have to be adjusted to each particular cell. This is done by partitioning each cell into triangles and by using [John Burkardt's implementation of the Dunavant rules](https://people.sc.fsu.edu/~jburkardt/cpp_src/triangle_dunavant_rule/triangle_dunavant_rule.html). The choice was also made not to pre-compute all quadrature rules for all cells and edges, but rather to compute them -- with a locally chosen degree of exactness -- when needed in the code. To reduce the computational cost, quadrature rules -- and the values of basis functions at quadrature nodes -- should only be computed once when looping over each cell, before being used, e.g., to form mass matrices.

The [HybridCore](@ref HArDCore2D::HybridCore) structure provides routines to do that. The method [cell_qrule(iT,doe)](@ref HArDCore2D::HybridCore#cell_qrule) calculates quadrature nodes and weights, exact up to the polynomial degree `doe` for an integration over cell number `iT`; see [edge_qrule(iE,doe)](@ref HArDCore2D::HybridCore#edge_qrule) for the equivalent over an edge. This quadrature rule, stored for example in `quadTE`, is then be provided to [basis_quad(type,i,quadTE,deg)](@ref HArDCore2D::HybridCore#basis_quad) which computes the values of the basis functions in cell/edge (depending on `type`=T or some other character) number `i` up to the specified degree `deg`; see also [grad_basis_quad(i,quadTE,deg)](@ref HArDCore2D::HybridCore#grad_basis_quad) to compute the gradients of the basis functions at the quadrature nodes.

Typically, these values are then passed on to [gram_matrix(f_quad,g_quad,Nf,Ng,quadTE,sym,weight)](@ref HArDCore2D::HybridCore#gram_matrix) which computes a ``generalised'' Gram matrix \f$(\int weight*f_i*g_j)_{ij}\f$ of the families of functions \f$(f_1,\ldots,f_{Nf})\f$ and \f$(g_1,\ldots,g_{Ng})\f$, provided at the quadrature notes by `f_quad` and `g_quad` (here, `sym` is a boolean indicating if the matrix is expected to be symmetric). This Gram matrix method is useful to compute mass and stiffness matrices.

Here is an example.

\code{.cpp}
// Create quadrature rule on cell `iT`. Here, `hho` is an instance of the `HybridCore` class. The degree of 
// exactness ensures that the rule will be exact for polynomial functions up to degree \f$K+L+1\f$
std::vector<HybridCore::qrule> quadT = hho.cell_qrule(iT, hho.Ldeg()+hho.K()+1);

// Compute values of basis functions, up to degree \f$(K+1)\$, at the quadrature nodes
std::vector<Eigen::VectorXd> phi_quadT = hho.basis_quad('T', iT, quadT, hho.nhighorder_dofs());

// Compute values of gradients of basis functions, as well as `diff` times these gradients, at the quadrature nodes
std::vector<Eigen::MatrixXd> dphiT_quadT = hho.grad_basis_quad(iT, quadT, hho.nhighorder_dofs());

// Create the mass matrix of basis functions up to degree \f$L\f$ and the stiffness matrix of the gradients
// kappa is a function that returns the diffusion tensor at the considered location
Eigen::MatrixXd MTT = hho.gram_matrix(phi_quadT, phi_quadT, hho.nlocal_cell_dofs(), hho.nhighorder_dofs(), quadT, true);
Eigen::MatrixXd StiffT = hho.gram_matrix(dphiT_quadT, dphiT_quadT, hho.nhighorder_dofs(), hho.nhighorder_dofs(), quadT, true, kappa);

// Grab the global index of the first edge of cell iT, compute quadrature nodes on this edge
size_t iF = mesh->cell(iT)->edge(0)->global_index();
std::vector<HybridCore::qrule> quadF = hho.edge_qrule(iF, 2*hho.K()+2);

// Compute the values of the cell basis function, and the edge basis function, at the quadrature points
// on the edge, and create the 'mass matrix' of cell-edge basis functions on the edge
std::vector<Eigen::VectorXd> phiT_quadF = hho.basis_quad('T', iT, quadF, hho.nhighorder_dofs());
std::vector<Eigen::VectorXd> phiF_quadF = hho.basis_quad('F', iF, quadF, hho.nlocal_edge_dofs());
Eigen::MatrixXd	MFT = hho.gram_matrix(phiF_quadF, phiT_quadF, hho.nlocal_edge_dofs(), hho.nhighorder_dofs(), quadF, false);
\endcode

At present, the quadrature rules available in the code support a total degree of exactness in the cells or on the edges up to 20.



<!--<a name="hmm">-->
<!--\section example Example - The Hybrid Mimetic Mixed (HMM) Scheme-->
<!--</a>-->

<!--We present in this section, an example implementation of a scheme. The following code assembles the matrices for the Hybrid Mimetic Mixed (HMM) scheme. The key part of the assembly is the computation of the local matrices \f$\mathbb{W}_K\f$ from the mesh data. As usual, we use the [Eigen](http://eigen.tuxfamily.org/) library for the matrices. The complete implementation can be found in the file `Examples/HMM/hmm.cpp`.-->

<!--\code{.cpp}-->
<!--// ---------------------------------------------------------------------------->
<!--//                        Create the HHO data structure-->
<!--// ---------------------------------------------------------------------------->

<!--// the mesh-->
<!--auto mesh_ptr = std::make_shared<EigenMesh3D>();-->
<!--mesh_3Dv& mesh = *mesh_ptr;-->

<!--// read from tetgen files-->
<!--ExtFileInput input_mesh ;-->
<!--input_mesh.TETGEN_format(mesh, mesh_file, offset);-->

<!--// Create the HHO-->
<!--HybridCore hho(mesh_ptr, 0, 0, 2, 2);-->

<!--// ---------------------------------------------------------------------------->
<!--//                          HMM parameters-->
<!--// ---------------------------------------------------------------------------->

<!--// Identity tensor-->
<!--Eigen::MatrixXd Lambda = Eigen::MatrixXd::Identity(3,3);-->

<!--// Exact solution-->
<!--auto u_exact = [](auto x, auto y, auto z) {-->
<!--  return sin(pi*x) * sin(pi*y) * sin(pi*z);-->
<!--};-->

<!--// Source term-->
<!--auto f = [](auto x, auto y, auto z) {-->
<!--  return 3.0 * pi * pi * sin(pi*x) * sin(pi*y) * sin(pi*z);-->
<!--};-->

<!--const double beta_K = 1;  // A "positive parameter"-->

<!--// ---------------------------------------------------------------------------->
<!--//                          Scheme Assembly-->
<!--// ---------------------------------------------------------------------------->

<!--const int ncells = mesh.n_region();   // Number of cells / 3D regions-->
<!--const int nfaces = mesh.n_face();     // Number of faces / 2D regions-->
<!--const int nbfaces = mesh.n_bface();   // Number of boundary faces-->

<!--const int ntotal_dofs = ncells + nfaces;-->

<!--// Matrices-->
<!--Eigen::SparseMatrix<double> A(ntotal_dofs, ntotal_dofs);-->
<!--Eigen::VectorXd b = Eigen::VectorXd::Zero(ntotal_dofs);-->
<!--std::vector<Eigen::Triplet<double>> triplets;-->

<!--// Fix the boundary faces to zero-->
<!--for (int i = 0; i < nbfaces; i++) {-->
<!--  const int iF = ncells + mesh.get_bnd_face(i);-->
<!--  triplets.emplace_back(iF, iF, 1);-->
<!--}-->

<!--// For each cell / region R in the mesh-->
<!--for (int iR = 0; iR < ncells; iR++) {-->
<!--  const int nlocalfaces = mesh.n_regn_face(iR);-->

<!--  // Cell barycenter and measure-->
<!--  auto xR = mesh.region_barycenter(iR);-->
<!--  auto mR = mesh.get_regn_measure(iR);-->

<!--  // Compute G_K, X_K, D_K-->
<!--  Eigen::MatrixXd G_K(d, nlocalfaces);-->
<!--  Eigen::MatrixXd X_K(nlocalfaces, d);-->
<!--  Eigen::MatrixXd D_K = Eigen::MatrixXd::Zero(nlocalfaces, nlocalfaces);-->

<!--  for (int ilF = 0; ilF < nlocalfaces; ilF++) {-->
<!--    const int iF = mesh.regn_face(iR, ilF);-->

<!--    // Face barycenter and measure-->
<!--    auto xF = mesh.face_barycenter(iF);-->
<!--    auto mF = mesh.get_face_measure(iF);-->

<!--    // Get the normal pointing out of iR from ilF-->
<!--    auto nRF = mesh.exterior_face_normal(iR, ilF);-->

<!--    // Fill G_K and X_K-->
<!--    G_K.col(ilF) = -mF / mR * nRF;-->
<!--    X_K.row(ilF) = (xF - xR).transpose();-->

<!--    // Fill in diagonal of D_K-->
<!--    auto dk_sigma = (xF - xR).dot(nRF);-->
<!--    D_K(ilF, ilF) = mF / dk_sigma * (K * nRF).dot(nRF);-->
<!--  }-->

<!--  // Construct W_K-->
<!--  Eigen::MatrixXd T_K = -Eigen::MatrixXd::Identity(nlocalfaces,nlocalfaces) - X_K * G_K;-->
<!--  Eigen::MatrixXd W_K = mR * G_K.transpose() * K * G_K + beta_K * beta_K * T_K.transpose() * D_K * T_K;-->

<!--  for (int ilF = 0; ilF < nlocalfaces; ilF++) {-->
<!--    const int iF = mesh.regn_face(iR, ilF);-->
<!--    for (int jlF = 0; jlF < nlocalfaces; jlF++) {-->
<!--      const int jF = mesh.regn_face(iR, jlF);-->

<!--      triplets.emplace_back(iR,iR, W_K(ilF, jlF));-->
<!--      triplets.emplace_back(iR, ncells + jF, -W_K(ilF, jlF));-->

<!--      if (mesh.is_internal_face(iF)) {-->
<!--        triplets.emplace_back(ncells + iF, iR, W_K(ilF, jlF));-->
<!--        triplets.emplace_back(ncells + iF, ncells + jF, -W_K(ilF, jlF));-->
<!--      }-->
<!--    }-->
<!--  }-->

<!--  // Assemble the source term-->
<!--  b(iR) = hho.integrate_over_region(iR, f);-->
<!--}-->

<!--// Assemble the sparse matrix A-->
<!--A.setFromTriplets(std::begin(triplets), std::end(triplets));-->
<!--\endcode-->

<!--An example implementation of a much more complicated scheme can be found in `Examples/HHO/DiffusionEquation.hpp`.-->


