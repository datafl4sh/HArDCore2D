// Core data structures and methods required to implement hybrid schemes in 2D (polynomial unknowns
// in the cells and on the edges, such as Hybrid High-order (HHO) schemes).
//
// Provides:
//  - Hybrid polynomial basis functions (on the cells and edges of the mesh)
//  - Generic routines to create quadrature nodes over cells and edges of the mesh
//  - Interpolation of general functions onto the HHO space
//  - Methods for integrating, evaluating, and computing norms of HHO solutions
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef HYBRIDCORE_HPP
#define HYBRIDCORE_HPP

#include <cassert>
#include <cmath>

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <mesh.hpp>
#include <cell.hpp>
#include <edge.hpp>

/*!	
* @defgroup HybridCore 
* @brief Classes providing tools to implement schemes having polynomial unknowns on the edges and in the cells
*/

namespace HArDCore2D {


/*!
*	\addtogroup HybridCore
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/** The HybridCore class provides convenient interfaces for performing
* integration over mesh cells and faces and handling polynomial basis
* functions

**/
/** The class also provides convenient interfaces for dealing with solutions to
   Hybrid High-Order schemes, such as the
    computation of integrals, norms and interpolants in the HHO space. */

class HybridCore {

public:
	///@brief Class constructor: initialises the data structure with the given mesh, and desired polynomial degrees of the basis functions
  HybridCore(
						const Mesh* mesh_ptr, ///< A pointer to the loaded mesh 
            const size_t K, ///< The degree of the edge polynomials 
						const size_t L ///< The degree of the cell polynomials 
             ); 

  // Basis functions
  using cell_basis_type = std::function<double(double, double)>;		///< type for cell basis
  using cell_gradient_type = std::function<Eigen::Vector2d(double, double)>;   ///< type for gradients of cell basis
	using edge_basis_type = std::function<double(double,double)>;   ///< type for edge basis
	using tensor_function_type = std::function<Eigen::Matrix2d(double,double)>;   ///< type for 2D tensors basis

	///@brief description of one node and one weight from a quadrature rule
  struct qrule {			
      double x, y, w;
      qrule(double x, double y, double w) : x(x), y(y), w(w) {}
  };

  /// Compute the size of the basis of 2-variate polynomials up to degree m
  size_t dim_Pcell(
										const size_t m /**< The maximum degree of the polynomial */
                   ) const;

	/// Compute the size of the basis of 1-variate polynomials up to degree m
  size_t dim_Pedge(
										const size_t m /**< The maximum degree of the polynomial */
                   ) const;

  /// Return a reference to the i'th basis function of the cell iT
  const cell_basis_type& cell_basis(
      size_t iT, /**< The global region number of the region */
      size_t i   /**< The index of the desired basis function */
      ) const;
	
	/// Return a reference to the i'th basis function of the edge iF
	const edge_basis_type& edge_basis(
      size_t iF, /**< The global edge number of the edge */
      size_t i   /**< The index of the desired basis function */
      ) const;
	
  /// Return a reference to the gradient of the i'th basis function of the celliT
  /** Note that the gradient functions are indexed the same as the basis
     functions. In particular, this
      means that the first gradient function will always be identically
     zero, as it is the gradient of the constant basis function. */
  const cell_gradient_type& cell_gradient(
      size_t iT, /**< The global region number of the region */
      size_t i   /**< The index of the desired basis function */
      ) const;
  const Mesh* get_mesh_ptr();		///< returns a pointer to the mesh

	/// Compute a quadrature rule of a given degree of exactness on a cell
  std::vector<HybridCore::qrule> cell_qrule(const size_t iT,		/**< Global index of the cell */
                                            const size_t doe		/**< Required degree of exactness */
																					) const;		/**< @returns list of quadrature nodes and weights */

	/// Compute a quadrature rule of a given degree of exactness on an edge
  std::vector<HybridCore::qrule> edge_qrule(const size_t iE,	/**< Global index of the edge */
                                            const size_t doe		/**< Required degree of exactness */
																					) const;  /**< @returns list of quadrature nodes and weights */

	Eigen::VectorXd restr(const Eigen::VectorXd &Xh, size_t iT) const;	///< Extract from a global vector Xh of unknowns the unknowns corresponding to cell iT
	double L2norm(const Eigen::VectorXd &Xh) const; ///< Compute L2 norm of a discrete function (using cell values)
	double H1norm(const Eigen::VectorXd &Xh) const; ///< Compute discrete H1 norm of a discrete function (using cell values)
	double Linf_edge(const Eigen::VectorXd &Xh) const;	///< Compute maximum of the coefficients on the edge basis functions

	/// Compute the interpolant in the discrete space of a continuous function
	template<typename Function>
	Eigen::VectorXd interpolate(const Function& f, size_t doe) const;

	/// Create the matrix of L2 products of two families (f_i) and (g_j) of functions
	/// (this is not really a Gram matrix, unless the two families are the same)
	Eigen::MatrixXd gram_matrix(
			const std::vector<Eigen::ArrayXd>& f_quad, ///< Values of functions (f1,f2,...) at the quadrature nodes 
			const std::vector<Eigen::ArrayXd>& g_quad, ///< Values of functions (g1,g2,...) at the quadrature nodes 
			const size_t& nrows, ///< Number of rows of the matrix - typically number of functions f_i (but could be less) 
			const size_t& ncols, ///< Number of columns of the matrix - typically number of functions g_j (but could be less)
			const std::vector<HybridCore::qrule>& quad, 	///< Quadrature nodes for integration 
			const bool& sym,		///< True if the matrix is pseudo-symmetric (that is, #f<=#g and f_i=g_i if i<=#f) 
			std::vector<double> L2weight = {} 	///< Optional weight for the L2 product. If provided, should be a std::vector<double> of the weight at the quadrature nodes
	) const;	/**< @returns The matrix \f$(\int f_i g_j)_{i=1\ldots nrows; j=1\ldots ncols}\f$ */

	/// Overloaded version of the previous one for vector-valued functions: the functions (F_i) and (G_j) are vector-valued functions
	Eigen::MatrixXd gram_matrix(
			const std::vector<Eigen::ArrayXXd>& F_quad,		///< Values of functions (F1,F2,...) at the quadrature nodes 
			const std::vector<Eigen::ArrayXXd>& G_quad,		///< Values of functions (G1,G2,...) at the quadrature nodes 
			const size_t& nrows,		///< Number of rows of the matrix - typically number of functions F_i (but could be less) 
			const size_t& ncols,		///< Number of rows of the matrix - typically number of functions G_j (but could be less) 
			const std::vector<HybridCore::qrule>& quad,		///< Quadrature nodes for integration 
			const bool& sym,		///< True if the matrix is pseudo-symmetric (that is, #F<=#G and F_i=G_i if i<=#F) 
			std::vector<Eigen::Matrix2d> L2Weight = {}	///< Optional weight for the L2 product. If provided, should be a std::vector<Eigen::Matrix2d> of the weight at the quadrature nodes
		) const ;  /**< @returns The matrix \f$(\int F_i \cdot G_j)_{i=1\ldots nrows; j=1\ldots ncols}\f$ */
	
	
	Eigen::VectorXd compute_weights(size_t iT) const; ///< Weights to compute cell unknowns from edge unknowns when l=-1
	
	/// Computes (cell or edge) basis functions at the quadrature nodes
	std::vector<Eigen::ArrayXd> basis_quad(
			const char B, 			///< T for cell, E for edge
			const size_t iTF, 		///<	global index of the cell/edge
			const std::vector<qrule> quad, 	///< quadrature nodes and weights on the cell/edge
			const size_t nbasis		///< number of basis functions to consider
	) const;	///< @returns phi_quad[i] = array listing values of phi_i at the quadrature nodes

	/// Compute \f$(\nabla \phi_i)_{i\in I}\f$ at the quadrature nodes, where \f$(\phi_i)_{i\in I}\f$ are the cell basis functions
	std::vector<Eigen::ArrayXXd> grad_basis_quad(
			const size_t iT, 											///< global index of the cell
			const std::vector<qrule> quad, 				///< quadrature rules in the cell
			const size_t nbasis										///< number of basis functions \f$phi_i\f$ to consider
	) const; ///< @returns dphi_quad[i]: array of size 2*nbq (where nbq=nb of quadrature nodes), with each column being \f$\nabla \phi_i\f$ at the corresponding quadrature node 


	/// Computes the integral of a discrete function using cell unknowns
	double integral(const Eigen::VectorXd &XTF) const; 
	/// Evaluates a discrete function in the cell iT at point (x,y)
	double evaluate_in_cell(const Eigen::VectorXd XTF, size_t iT, double x, double y) const; 
	/// Evaluates a discrete function on the edge iF at point (x,y)
	double evaluate_in_edge(const Eigen::VectorXd XTF, size_t iF, double x, double y) const;

	const size_t K();		///< polynomial degree of edge unknowns
	const int L();		///< polynomial degree of cell unknowns
	const int Ldeg();		///< usually equal to L, but put at 0 if L=-1
  inline size_t ntotal_dofs();	///< Total number of degrees of freedom
	inline size_t nlocal_cell_dofs(); ///< number of degrees of freedom in each cell (dimension of polynomial space)
  inline size_t ntotal_cell_dofs(); ///< total number of cell degrees of freedom
	inline size_t nlocal_edge_dofs();	///< number of degrees of freedom on each edge (dimension of polynomial space)
	inline size_t ntotal_edge_dofs();	///< total number of edge degrees of freedom
	inline size_t ninternal_edge_dofs();	///< total number of edge degrees of freedom for internal edges
	inline size_t nboundary_edge_dofs();	///< total number of edge degrees of freedom for boundary edges
	inline size_t nhighorder_dofs();	///< total number of cell degrees of freedom with polynomials up to order k+1
	inline size_t ngradient_dofs();  ///< total number of degrees of freedom for gradients

	/// To integrate a function over a cell
  template <typename Function>
  void quadrature_over_cell(const size_t iT, const Function& f) const;
	/// To integrate a function over an edge
	template <typename Function>
	void quadrature_over_edge(const size_t iF, const Function& f) const;
	///Integrates a function over a cell. Use with parcimony, expensive (re-compute quadratures)
	template <typename Function>
	double integrate_over_cell(const size_t iT, const Function& f) const;
	///Integrates a function over an edge. Use with parcimony, expensive (re-compute quadratures)
	template <typename Function>
	double integrate_over_edge(const size_t iF, const Function& f) const;
	///Integrates a function over the domaine. Use with parcimony, expensive (re-compute quadratures)
	template <typename Function>
	double integrate_over_domain(const Function& f) const;

private:
	std::pair<std::vector<HybridCore::cell_basis_type>,
		        std::vector<HybridCore::cell_gradient_type> >
  create_cell_basis(size_t iG) const;			///< creates_the basis for cell iG of degree K.
	///@brief creates the basis for edge iF of degree L.
	std::vector<HybridCore::edge_basis_type> create_edge_basis(size_t iF) const;
	// Mesh data
	const Mesh* _mesh_ptr;  // Shared ptr to mesh data
	// const Mesh& _mesh;  // Const-reference to the mesh data for
	//                      // immutable access
	/// Degree of the edge polynomials
	const size_t _K;
	/// Degree of the cell polynomials (can be -1).
	const int _L;
	const size_t _Ldeg;
	/// Number of degrees of freedom per cell polynomial
	const size_t _nlocal_cell_dofs;
	/// Number of degrees of freedom per cell for high-order (K+1) polynomials
	//const int _nhighorder_dofs;
	/** Number of degrees of freedom in the space of gradients of high-order
	 cell polynomials*/

	const size_t _nlocal_edge_dofs;
	/// Number of degrees of freedom per cell for high-order (K+1) polynomials
	const size_t _nhighorder_dofs;
	/// Number of degrees of freedom in the space of gradients of high-order cell polynomials
	const size_t _ngradient_dofs;
	/// Total number of cell degrees of freedom over the entire mesh
	const size_t _ntotal_cell_dofs;
	/// Total number of edge degrees of freedom over the entire mesh
	const size_t _ntotal_edge_dofs;
	/// Total number of edge degrees of freedom inside the domain
	const size_t _ninternal_edge_dofs;
	/// Total number of edge degrees of freedom on the boundary of the domain
	const size_t _nboundary_edge_dofs;

	/// Total number of degrees of freedom over the entire mesh
	const size_t _ntotal_dofs;

	// Collections of basis functions for each cell and face
	std::vector<std::vector<cell_basis_type> > _cell_bases;
	std::vector<std::vector<cell_gradient_type> > _cell_gradients;
	std::vector<std::vector<edge_basis_type> > _edge_bases;

	// offset for quadrature rules, should be 0 except for testing purposes
	const int _offset_doe = 0;	

};



// ----------------------------------------------------------------------------
//                            Implementation
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ------- 'Easy' integration routines, expensive (re-compute quad rules), use as little as possible


template <typename Function>
void HybridCore::quadrature_over_cell(const size_t iT, const Function& f) const {
  assert(iT < _mesh_ptr->n_cells());
  std::vector<HybridCore::qrule> quadT = cell_qrule(iT, 2 * _Ldeg + 5);
  for (size_t iqn = 0; iqn < quadT.size(); iqn++) {
      f(iqn, quadT[iqn].x, quadT[iqn].y, quadT[iqn].w);
  }
}

template <typename Function>
void HybridCore::quadrature_over_edge(const size_t iF, const Function& f) const {
  assert(iF < _mesh_ptr->n_edges());

  std::vector<HybridCore::qrule> quadF = edge_qrule(iF, 2 * _K + 5);
  for (size_t iqn = 0; iqn < quadF.size(); iqn++) {
     f(iqn, quadF[iqn].x, quadF[iqn].y, quadF[iqn].w);
  }
}

template <typename Function>
double HybridCore::integrate_over_cell(const size_t iT, const Function& f) const {
  assert(iT < _mesh_ptr->n_cells());
  double ans = 0.0;
  quadrature_over_cell(iT, [&ans, &f](size_t iQN, double x, double y,
                                      double w) { ans += w * f(x, y); });
  return ans;
}

template <typename Function>
double HybridCore::integrate_over_edge(const size_t iF, const Function& f) const {
  assert(iF < _mesh_ptr->n_edges());
  double ans = 0.0;
  quadrature_over_edge(iF, [&ans, &f](size_t iQN, double x, double y,
                                      double w) { ans += w * f(x, y); });
  return ans;
}

template <typename Function>
double HybridCore::integrate_over_domain(const Function& f) const{
	double value = 0.0;
	for (size_t iT = 0; iT < _mesh_ptr->n_cells(); iT++){
		value += integrate_over_cell(iT, f);
	}
	return value;
}





// ----------------------------------------------------------
// ------- Interpolates continuous function on discrete space

template<typename Function>
Eigen::VectorXd HybridCore::interpolate(const Function& f, size_t doe) const {
  Eigen::VectorXd XTF = Eigen::VectorXd::Zero(_ntotal_dofs);

  for (size_t iT = 0; iT < _mesh_ptr->n_cells(); iT++) {

    // Mass matrix in cell
		std::vector<HybridCore::qrule> quadT = cell_qrule(iT, doe);
		std::vector<Eigen::ArrayXd> phi_quadT = basis_quad('T', iT, quadT, _nlocal_cell_dofs);
		Eigen::MatrixXd MT = gram_matrix(phi_quadT, phi_quadT, _nlocal_cell_dofs, _nlocal_cell_dofs, quadT, true);

		// Vector of integrals of f against basis functions
    Eigen::VectorXd bT = Eigen::VectorXd::Zero(_nlocal_cell_dofs);
    for (size_t i = 0; i < _nlocal_cell_dofs; i++) {
      const auto& phi_i = cell_basis(iT, i);
      bT(i) += integrate_over_cell(iT, [&phi_i, &f](double x, double y) {
        return phi_i(x,y) * f(x,y);
      });
    }

		// Vector of coefficients (on cell basis functions) of the L2(T) projection of f
    Eigen::VectorXd UT = MT.ldlt().solve(bT);

		// Fill in the complete vector of DOFs
    size_t offset_T = iT * _nlocal_cell_dofs;
    XTF.segment(offset_T, _nlocal_cell_dofs) = UT;

    for (size_t ilF = 0; ilF < _mesh_ptr->cell(iT)->n_edges(); ilF++) {
      size_t iF = _mesh_ptr->cell(iT)->edge(ilF)->global_index();

			// Mass matrix on face
			std::vector<HybridCore::qrule> quadF = edge_qrule(iF, doe);
			std::vector<Eigen::ArrayXd> phiF_quadF = basis_quad('F', iF, quadF, _nlocal_edge_dofs);
			Eigen::MatrixXd MF = gram_matrix(phiF_quadF, phiF_quadF, _nlocal_edge_dofs, _nlocal_edge_dofs, quadF, true);

			// Vector of integral of f against face basis functions
	    Eigen::VectorXd bF = Eigen::VectorXd::Zero(_nlocal_edge_dofs);
      for (size_t i = 0; i < _nlocal_edge_dofs; i++) {
        const auto& phi_i = edge_basis(iF, i);
        bF(i) += integrate_over_edge(iF, [&phi_i, &f](double x, double y) {
          return phi_i(x,y) * f(x,y);
        });
      }

			// Vector of coefficients (on cell basis functions) of the L2(F) projection of f
		  Eigen::VectorXd UF = MF.ldlt().solve(bF);

			// Fill in he complete vector of DOFs
	    XTF.segment(_ntotal_cell_dofs + iF * _nlocal_edge_dofs, _nlocal_edge_dofs) = UF;
    }

		// Special case of L=-1, we replace the previously computed cell value with the proper average of face values
		if ( _L==-1 ){
			Eigen::VectorXd barycoefT = compute_weights(iT);
			XTF(iT) = 0;
	    for (size_t ilF = 0; ilF < _mesh_ptr->cell(iT)->n_edges(); ilF++) {
				size_t iF = _mesh_ptr->cell(iT)->edge(ilF)->global_index();
				XTF(iT) += barycoefT(ilF) * XTF(_ntotal_cell_dofs + iF);
			}
		}

	}

  return XTF;
}



// --------------------------------------------------------------------------------------------------
// ------- Functions that return class elements

size_t HybridCore::nlocal_cell_dofs() { return _nlocal_cell_dofs; }
size_t HybridCore::ntotal_cell_dofs() { return _ntotal_cell_dofs; }
size_t HybridCore::nlocal_edge_dofs() { return _nlocal_edge_dofs; }
size_t HybridCore::ntotal_edge_dofs() { return _ntotal_edge_dofs; }
size_t HybridCore::nhighorder_dofs() { return _nhighorder_dofs; }
size_t HybridCore::ngradient_dofs() { return _ngradient_dofs; }

size_t HybridCore::ninternal_edge_dofs() { return _ninternal_edge_dofs; }
size_t HybridCore::nboundary_edge_dofs() { return _nboundary_edge_dofs; }
size_t HybridCore::ntotal_dofs() { return _ntotal_dofs; }

//@}

}  // end of namespace HArDCore2D

#endif /* HYBRIDCORE_HPP */
