// Implementation of the HHO scheme in 2D for the diffusion equation
//
//   { -div(K \grad(u)) = f,       inside Omega
//   { K \grad(u) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
// 
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef _DIFFUSIONEQUATION_HPP
#define _DIFFUSIONEQUATION_HPP

#include <functional>
#include <utility>
#include <iostream>

#include <boost/timer.hpp>
#include <stdlib.h>     /* exit, EXIT_FAILURE */

// Matrices and linear solvers
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include "Eigen/MA41.cpp"


#include "mesh2d.hpp"
#include "hybridcore2d.hpp"
#include "quad2d.hpp"

/*!
* @defgroup HHO_diffusion
* @brief HHO scheme for diffusion equation -div(Diff grad u)=f
* @brief The scheme here corresponds to a variable diffusion that is piecewise constant in the cells
*/

namespace HArDCore2D {

/*!
*	@addtogroup HHO_diffusion
* @{
*/
// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The DiffusionEquation class provides tools to implement the HHO method for the diffusion problem

class DiffusionEquation {

// Types
public:
  using scalar_function_type = std::function<double(double,double)>;		///< type for function R^2->R
  using vector_function_type = std::function<Eigen::VectorXd(double,double)>;		///< type for function R^2->R^2
  using tensor_function_type = std::function<Eigen::Matrix2d(double,double)>;		///< type for function R^2->R^{2x2}

	///@brief Constructor of the class
  DiffusionEquation(
		tensor_function_type kappa, 	///< diffusion tensor
		scalar_function_type source,  ///< source term
		size_t BC, 										///< type of boundary conditions (0 for Dirichlet, 1 for Neumann)
		scalar_function_type exact_solution, 	///< exact solution
		vector_function_type grad_exact_solution, 	///< gradient of the exact solution
		std::string solver_type		///< type of solver to use for the global system (bicgstab at the moment)
		);

  /// Assemble and solve the scheme
  Eigen::VectorXd solve(HybridCore& hho);

	/// Discrete H1 norm of an hybrid function
  double H1Norm(HybridCore& hho, const Eigen::VectorXd XTF); 

  /// From a hybrid function, computes a vector of values at the vertices of the mesh
  Eigen::VectorXd VertexValues(
			HybridCore& hho, 	///< current instance of the hybridcore class
			const Eigen::VectorXd XTF, 		///< hybrid function (cell and edge polynomials)
			const char from_dofs		///< Type of unknowns to use: cell (T) or edge (E)
			);

  double get_assembly_time() const; ///< cpu time to assemble the scheme
  double get_solving_time() const;	///< cpu time to solve the scheme
  double get_solving_error() const;	///< residual after solving the scheme
  double get_itime(size_t idx) const;		///< various intermediate assembly times

private:
  /// Compute the local diffusion operator in the cell iT
  Eigen::MatrixXd diffusion_operator(HybridCore& hho, const size_t iT) const;

  /// Compute the local load operator (the source term) in the cell iT
  Eigen::VectorXd load_operator(HybridCore& hho, const size_t iT) const;

  const tensor_function_type kappa;
  const scalar_function_type source;
	const size_t BC;
  const scalar_function_type exact_solution;
  const vector_function_type grad_exact_solution;
	const std::string solver_type;

  // Computation statistics
  double assembly_time;
  double solving_time;
  double solving_error;
	mutable std::vector<double> itime = std::vector<double>(5, 0.0);

};

DiffusionEquation::DiffusionEquation(tensor_function_type kappa, scalar_function_type source, size_t BC, scalar_function_type exact_solution, vector_function_type grad_exact_solution, std::string solver_type)
  : kappa(std::move(kappa)),
    source(std::move(source)),
		BC(std::move(BC)),
		exact_solution(std::move(exact_solution)),
		grad_exact_solution(std::move(grad_exact_solution)),
		solver_type(std::move(solver_type)) {
  // Do nothing
}

Eigen::VectorXd DiffusionEquation::solve(HybridCore &hho) {

	boost::timer timer;  // Time the matrix assembly
  const auto mesh = hho.get_mesh_ptr();

	// System matrix
  Eigen::SparseMatrix<double> A(hho.ntotal_edge_dofs(), hho.ntotal_edge_dofs());
	// If static condensation (L>=0): matrix to recover cell unknowns
	// If barycentric elimination (L=-1): matrix to recover cell unknowns
  Eigen::SparseMatrix<double> A_sc_be(hho.ntotal_cell_dofs(), hho.ntotal_edge_dofs());

  Eigen::VectorXd B = Eigen::VectorXd::Zero(hho.ntotal_edge_dofs());
  Eigen::VectorXd B_sc = Eigen::VectorXd::Zero(hho.ntotal_cell_dofs());

  //------------------------------------------------------------------------------
  // Assemble local contributions

  // Global quadrature rule for the cells

  Eigen::VectorXd cell_quadrature = Eigen::VectorXd::Zero(hho.ntotal_cell_dofs());


  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<Eigen::Triplet<double>> triplets_sc_be;    // Triplets to create A_sc or A_be

  auto total_measure = 0.0;

  for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
		Cell2D* iCell = mesh->cell(iT);
	
    total_measure += iCell->measure();

    // Total number of face degrees of freedom local to this cell (adjacent faces to the cell)
    size_t nlocal_edges = iCell->n_edges();
    size_t edge_dofs = nlocal_edges * hho.nlocal_edge_dofs();
    // Compute local operators

		Eigen::VectorXd bT = load_operator(hho, iT);
		Eigen::MatrixXd aT = diffusion_operator(hho, iT);

		// Local matrix and right-hand side on the face unknowns
		Eigen::MatrixXd MatF = Eigen::MatrixXd::Zero(edge_dofs,edge_dofs);
		Eigen::VectorXd bF = Eigen::VectorXd::Zero(edge_dofs);

		if (hho.L()>=0) {
			// STATIC CONDENSATION OF ELEMENT UNKNOWNS

		  // Perform static condensation
		  Eigen::MatrixXd ATT = aT.topLeftCorner(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs());
		  Eigen::MatrixXd ATF = aT.topRightCorner(hho.nlocal_cell_dofs(), edge_dofs);
		  Eigen::MatrixXd AFF = aT.bottomRightCorner(edge_dofs, edge_dofs);

		  Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
		  invATT.compute(ATT);
		  
		   Eigen::MatrixXd invATT_ATF = invATT.solve(ATF);
		  Eigen::VectorXd invATT_bTcell = invATT.solve(bT.head(hho.nlocal_cell_dofs()));
		  MatF = AFF - ATF.transpose() * invATT_ATF;
		  
		  bF = bT.tail(edge_dofs) - ATF.transpose() * invATT_bTcell;
						
		  // Assemble static condensation operator
		  B_sc.segment(iT * hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs()) = invATT_bTcell;
		  for (size_t i = 0; i < hho.nlocal_cell_dofs(); i++) {
		    for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
		      const size_t jF = iCell->edge(jlF)->global_index();
		      for (size_t jk = 0; jk < hho.nlocal_edge_dofs(); jk++) {
		        const size_t jLocal = jlF * hho.nlocal_edge_dofs() + jk;
		        const size_t jGlobal = jF * hho.nlocal_edge_dofs() + jk;
		        triplets_sc_be.emplace_back(iT * hho.nlocal_cell_dofs() + i, jGlobal, invATT_ATF(i, jLocal));
		      }
		    }
		  }
		} else {
			// BARYCENTRIC ELIMINATION OF ELEMENT UNKNOWNS
			// Create reduction matrix: 1+nlocal_edges * nlocal_edges matrix with the coefficients on the first row, and the identity below. When multiplied by the face unknowns, return cell and face unknowns
			Eigen::MatrixXd red_matT = Eigen::MatrixXd::Zero(1+nlocal_edges,nlocal_edges);
			red_matT.row(0) = hho.compute_weights(iT);
			red_matT.bottomRightCorner(nlocal_edges,nlocal_edges) = Eigen::MatrixXd::Identity(nlocal_edges,nlocal_edges);

			bF = red_matT.transpose() * bT;
			MatF = red_matT.transpose() * aT * red_matT;

			// To recover cell unknown
			for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
	      const size_t jF = iCell->edge(jlF)->global_index();
        const size_t jGlobal = jF * hho.nlocal_edge_dofs();
		  	triplets_sc_be.emplace_back(iT, jGlobal, red_matT(0,jlF));
			}

		}


    // Assemble local contribution into global matrix
    for (size_t ilF = 0; ilF < nlocal_edges; ilF++) {
      const size_t iF = iCell->edge(ilF)->global_index();
      for (size_t ik = 0; ik < hho.nlocal_edge_dofs(); ik++) {
        const size_t iLocal = ilF * hho.nlocal_edge_dofs() + ik;
        const size_t iGlobal = iF * hho.nlocal_edge_dofs() + ik;
        for (size_t jlF = 0; jlF < nlocal_edges; jlF++) {
          const size_t jF = iCell->edge(jlF)->global_index();
          for (size_t jk = 0; jk < hho.nlocal_edge_dofs(); jk++) {
            const size_t jLocal = jlF * hho.nlocal_edge_dofs() + jk;
            const size_t jGlobal = jF * hho.nlocal_edge_dofs() + jk;
            triplets.emplace_back(iGlobal, jGlobal, MatF(iLocal, jLocal));
          }
        }
        B(iGlobal) += bF(iLocal);
      }
    }
	
    // Record cell quadrature
		if (BC==1){
	    for (size_t i = 0; i < hho.nlocal_cell_dofs(); i++) {
  	    const auto& phi_i = hho.cell_basis(iT, i);
  	    cell_quadrature(iT * hho.nlocal_cell_dofs() + i) = hho.integrate_over_cell(iT, phi_i);
  	  }
		}

  }
	if (BC==0){
		// Dirichlet boundary conditions: the equations for boundary faces are replaced
		//	with MF UF = (int of phi_i u)_i. This forces UF to represent the coefficients
		// 	of the L2(F) projection of u on the basis (phi_i)_i
		std::vector<Edge2D *> b_edges = mesh->get_b_edges(); // List of boundary edges
		for (size_t ibF = 0; ibF < mesh->n_b_edges(); ibF++){
			size_t iF = b_edges[ibF]->global_index();
			// Mass matrix and boundary values
			std::vector<HybridCore::qrule> quadF = hho.edge_qrule(iF, 2*hho.K());
			std::vector<Eigen::VectorXd> phiF_quadF = hho.basis_quad('F', iF, quadF, hho.nlocal_edge_dofs());
			Eigen::MatrixXd MF = hho.gram_matrix(phiF_quadF, phiF_quadF, hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), quadF, true);
			for (size_t j = 0; j < hho.nlocal_edge_dofs(); j++){
				// Erase row in matrix and replace
				size_t jGlobal = iF * hho.nlocal_edge_dofs() + j;
				triplets.erase(std::remove_if(std::begin(triplets), std::end(triplets),
		              [&jGlobal](const auto& x) { return (x.row() == (int) jGlobal); }), std::end(triplets));
				for (size_t k = 0; k < hho.nlocal_edge_dofs(); k++){
					size_t kGlobal = iF * hho.nlocal_edge_dofs() + k;
					triplets.emplace_back(jGlobal, kGlobal, MF(j,k));
				}
				// Replace boundary value
        const auto& phi_j = hho.edge_basis(iF, j);
				B(jGlobal) = hho.integrate_over_edge(iF, [this, &phi_j](double x, double y) {
          return phi_j(x,y) * exact_solution(x,y); });
			}

		}
	} else if (BC==1){
		// Remove a row in the matrix and fix the first degree of freedom
		triplets.erase(std::remove_if(std::begin(triplets), std::end(triplets),
		                              [](const auto& x) { return (x.row() == 0); }), std::end(triplets));
		triplets.emplace_back(0, 0, 1);
		B(0) = 0;
	}
  // Assemble the linear system
  A.setFromTriplets(std::begin(triplets), std::end(triplets));
  A_sc_be.setFromTriplets(std::begin(triplets_sc_be), std::end(triplets_sc_be));

  assembly_time = timer.elapsed();  // Record final assembly time
  timer.restart();  // Time the solving

  // Solve the condensed system
	Eigen::VectorXd xF = Eigen::VectorXd::Zero(hho.ntotal_edge_dofs());
//	if (solver_type == "ma41") {
//		Eigen::MA41<Eigen::SparseMatrix<double>, Eigen::VectorXd> solver;
//		solver.analyzePattern(A);
//		solver.factorize(A);
//		xF = solver.solve(B);
//	} else {
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		xF = solver.solve(B);
		std::cout << "#iterations:     " << solver.iterations() << std::endl;
		std::cout << "estimated error: " << solver.error()      << std::endl;
//	}
  // Recover the cell unknowns
  Eigen::VectorXd XTFh = Eigen::VectorXd::Zero(hho.ntotal_dofs());
  XTFh.tail(hho.ntotal_edge_dofs()) = xF;

	if (hho.L()>=0) {
	  XTFh.head(hho.ntotal_cell_dofs()) = B_sc - A_sc_be * xF.head(hho.ntotal_edge_dofs());
	} else {
	  XTFh.head(hho.ntotal_cell_dofs()) = A_sc_be * xF.head(hho.ntotal_edge_dofs());
	}

	// Only Neumann: translate to get the proper average
	if (BC==1){
		// Compute average to translate
		auto average = cell_quadrature.dot(XTFh.head(hho.ntotal_cell_dofs())) / total_measure;
		auto average_exact_sol = hho.integrate_over_domain([&](auto x, auto y) {
			return exact_solution(x,y);
			}) / total_measure;

		// Translate the cells
		for (size_t iT = 0; iT < mesh->n_cells(); iT++) {
		  XTFh(iT * hho.nlocal_cell_dofs()) += average_exact_sol - average;
		}

		for (size_t iF = 0; iF < mesh->n_edges(); iF++) {
		  XTFh(hho.ntotal_cell_dofs() + iF * hho.nlocal_edge_dofs()) += average_exact_sol - average;
		}
	}

  solving_time = timer.elapsed();  // Record the final solving time

  // Compute the error in the linear system
  solving_error = (A * xF - B).norm();

  return XTFh;
}

//******************************** 
//		local diffusion matrix 
//********************************

Eigen::MatrixXd DiffusionEquation::diffusion_operator(HybridCore &hho, const size_t iT) const {

	boost::timer timeint;

  const auto mesh = hho.get_mesh_ptr();

  // Total number of degrees of freedom local to this cell (cell and its adjacent faces)
  size_t local_dofs = hho.nlocal_cell_dofs() + mesh->cell(iT)->n_edges() * hho.nlocal_edge_dofs();

  //-------------------  Initialisatons ------------------------------//

//timeint.restart();
	
	// Diffusion in the cell.
	tensor_function_type kappaT = [&](double x, double y){
		// Constant in the cell
		Eigen::Vector2d cellcenter = mesh->cell(iT)->center_mass();
			return kappa(cellcenter.x(), cellcenter.y());
		// Variable in the cell - the scheme may not provide optimal convergence rate if the diffusion is actually variable in the cell!
		// If the diffusion is piecewise constant, choosing the previous version ensures that, for edge integrals, it's the value inside the cell that is computed (essential in case kappa is discontinuous across the edges)
//			return kappa(x, y);
	};

	// Compute cell quadrature points and values of cell basis functions and gradients at these points
	std::vector<HybridCore::qrule> quadT = hho.cell_qrule(iT, hho.Ldeg()+hho.K()+1);
	std::vector<Eigen::VectorXd> phi_quadT = hho.basis_quad('T', iT, quadT, hho.nhighorder_dofs());
	std::vector<Eigen::MatrixXd> dphiT_quadT = hho.grad_basis_quad(iT, quadT, hho.nhighorder_dofs());

//itime[0] += timeint.elapsed();
//timeint.restart();

	// Preparing matrices to store face mass matrices (for face and cell basis functions), and face-cell mass matrices
	// (will be useful for the stabilisation term)
  std::vector<Eigen::MatrixXd> MFF(mesh->cell(iT)->n_edges(), Eigen::MatrixXd::Zero(hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs()));
  std::vector<Eigen::MatrixXd> MFT(mesh->cell(iT)->n_edges(), Eigen::MatrixXd::Zero(hho.nlocal_edge_dofs(), hho.nhighorder_dofs()));

//itime[1] += timeint.elapsed();
//timeint.restart();


  //-------------------- Compute PT, matrix of potential reconstruction ---------//

	// Stiffness matrix (kappaT dphi_i,dphi_j)_T for phi_i, phi_j up to degree K+1
	Eigen::MatrixXd StiffT = hho.gram_matrix(dphiT_quadT, dphiT_quadT, hho.nhighorder_dofs(), hho.nhighorder_dofs(), quadT, true, kappaT);  

	// Mass matrix of (phi_i,phi_j)_T for phi_i up to degree L and phi_j up to degree K+1
	Eigen::MatrixXd MTT = hho.gram_matrix(phi_quadT, phi_quadT, hho.nlocal_cell_dofs(), hho.nhighorder_dofs(), quadT, true);

	// Row vector of (phi_j,1)_T for phi_j up to degree K+1, and LT^t*LT (used later to impose average condition on PT)
	Eigen::VectorXd LT = (MTT.row(0)).transpose();
	Eigen::MatrixXd LTtLT = LT * (LT.transpose());

	// Right-hand side: we start with volumetric terms (dphi_i,dphi_j)_T for phi_i up to degree K+1 and phi_j up to degree L
  Eigen::MatrixXd BP = Eigen::MatrixXd::Zero(hho.nhighorder_dofs(), local_dofs);
	BP.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs()) = StiffT.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs());

  // Boundary terms in BP
  for (size_t ilF = 0; ilF < mesh->cell(iT)->n_edges(); ilF++) {
    // Offset for face unknowns
    const size_t offset_F = hho.nlocal_cell_dofs() + ilF * hho.nlocal_edge_dofs();
    const size_t iF = mesh->cell(iT)->edge(ilF)->global_index();
    const auto& nTF = mesh->cell(iT)->edge_normal(ilF);

		// Compute face quadrature points and values of basis functions (and gradients) at these points
		std::vector<HybridCore::qrule> quadF = hho.edge_qrule(iF, 2*hho.K()+2);
		size_t nbqF = quadF.size();
		
		std::vector<Eigen::VectorXd> phiT_quadF = hho.basis_quad('T', iT, quadF, hho.nhighorder_dofs());
		std::vector<Eigen::VectorXd> phiF_quadF = hho.basis_quad('F', iF, quadF, hho.nlocal_edge_dofs());
		std::vector<Eigen::MatrixXd> dphiT_quadF = hho.grad_basis_quad(iT, quadF, hho.nhighorder_dofs());

		// kappaT*nTF on quadrature points
		std::vector<Eigen::Vector2d> kappaT_nTF_quadF(nbqF, Eigen::Vector2d::Zero());
		for (size_t iqn = 0; iqn < nbqF; iqn++){
			kappaT_nTF_quadF[iqn] = kappaT(quadF[iqn].x, quadF[iqn].y)*nTF;
		}

    for (size_t i = 1; i < hho.nhighorder_dofs(); i++) {
			// We do not need i=0 because it corresponds to dphi_i=0
      for (size_t j = 0; j < hho.nlocal_cell_dofs(); j++) {
				for (size_t iqn = 0; iqn < nbqF; iqn++){
					BP(i,j) -= quadF[iqn].w * (dphiT_quadF[i].col(iqn)).dot(kappaT_nTF_quadF[iqn]) * phiT_quadF[j](iqn);
				}
      }
      for (size_t j = 0; j < hho.nlocal_edge_dofs(); j++) {
				for (size_t iqn=0; iqn < nbqF; iqn++){
					BP(i, offset_F + j) += quadF[iqn].w * (dphiT_quadF[i].col(iqn)).dot(kappaT_nTF_quadF[iqn]) * phiF_quadF[j](iqn);
				}
      }
    }

		// Face mass matrix, and cell-face mass matrix, for stabilisation term below
		MFF[ilF] = hho.gram_matrix(phiF_quadF, phiF_quadF, hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), quadF, true);
		MFT[ilF] = hho.gram_matrix(phiF_quadF, phiT_quadF, hho.nlocal_edge_dofs(), hho.nhighorder_dofs(), quadF, false);
  }

	// Compute PT
	double scalT = StiffT.trace() / std::pow(LT.norm(), 2);
	BP.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs()) += 
					scalT * LTtLT.topLeftCorner(hho.nhighorder_dofs(), hho.nlocal_cell_dofs());
	Eigen::MatrixXd PT = ((StiffT+scalT*LTtLT).ldlt()).solve(BP);

	// Consistent component (K \nabla pT, \nabla pT)_T in local bilinear form
  Eigen::MatrixXd ATF = PT.transpose() * StiffT * PT;

  //-------------------- Compute stabilisation term sT ---------//

  Eigen::MatrixXd STF = Eigen::MatrixXd::Zero(local_dofs, local_dofs);

	// Cell residual delta_T^l = pi_T^l (rT uT) - u_T
	Eigen::MatrixXd MTTKinv = ( MTT.topLeftCorner(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs()) ).inverse();
  Eigen::MatrixXd deltaTL = MTTKinv * MTT * PT;
  deltaTL.topLeftCorner(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs()) -=
        Eigen::MatrixXd::Identity(hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs());

  for (size_t ilF = 0; ilF < mesh->cell(iT)->n_edges(); ilF++) {
		auto hF = mesh->cell(iT)->edge(ilF)->measure();
    auto xF = mesh->cell(iT)->edge(ilF)->center_mass();

    auto kappa_TF = kappaT(xF.x(), xF.y()).trace();

		// Face residual delta_TF^k = pi_F^k (rT uT) - u_F
		Eigen::MatrixXd MFFinv = MFF[ilF].inverse();
    Eigen::MatrixXd deltaTFK = MFFinv * MFT[ilF] * PT;
    deltaTFK.block(0, hho.nlocal_cell_dofs() + ilF * hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs()) -=
        Eigen::MatrixXd::Identity(hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs());

		// Stabilisation term
		Eigen::MatrixXd dTFKminusdTL = deltaTFK - MFFinv * MFT[ilF].topLeftCorner(hho.nlocal_edge_dofs(), hho.nlocal_cell_dofs()) * deltaTL;

		STF += (kappa_TF / hF) * dTFKminusdTL.transpose() * MFF[ilF] *  dTFKminusdTL;

  }
//itime[2] += timeint.elapsed();
//timeint.restart();

	// Adjust local bilinear form with stabilisation term
  ATF += STF;

  return ATF;

}


//******************************** 
//		local load term 
//********************************

Eigen::VectorXd DiffusionEquation::load_operator(HybridCore &hho, const size_t iT) const {
	// Load for the cell DOFs (first indices) and face DOFs (last indices)
  const auto mesh = hho.get_mesh_ptr();
	size_t cell_edge_dofs = hho.nlocal_cell_dofs() + mesh->cell(iT)->n_edges()*hho.nlocal_edge_dofs();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(cell_edge_dofs);

	// Quadrature points and values of cell basis functions at these points
	std::vector<HybridCore::qrule> quadT = hho.cell_qrule(iT, 2*hho.K()+10);
	size_t nbq = quadT.size();
	std::vector<Eigen::VectorXd> phi_quadT = hho.basis_quad('T', iT, quadT, hho.nlocal_cell_dofs());

	// Value of source times quadrature weights at the quadrature points
	Eigen::VectorXd weight_source_quad = Eigen::VectorXd::Zero(nbq);
	for (size_t iqn = 0; iqn < nbq; iqn++){
		weight_source_quad(iqn) = quadT[iqn].w * source(quadT[iqn].x, quadT[iqn].y);
	}

	for (size_t i=0; i < hho.nlocal_cell_dofs(); i++){
		b(i) = weight_source_quad.dot(phi_quadT[i]);
	}
	// Boundary values, if we have a boundary cell
	if (mesh->cell(iT)->is_boundary()){
		if (BC==0){
			// Dirichlet BCs: no source terms on these edges
		} else if (BC==1) {
			// Neumann BCs
			for (size_t ilF = 0; ilF < mesh->cell(iT)->n_edges(); ilF++) {
			  const size_t iF = mesh->cell(iT)->edge(ilF)->global_index(); 
				// BC on boundary faces
				if (mesh->cell(iT)->edge(ilF)->is_boundary()){
				  // Offset for face unknowns
				  const size_t offset_F = hho.nlocal_cell_dofs() + ilF * hho.nlocal_edge_dofs();
					// Normal to the face
				  const auto& nTF = mesh->cell(iT)->edge_normal(ilF);
					// for each DOF of the boundary face
					for (size_t i = 0; i < hho.nlocal_edge_dofs(); i++){
				    const auto& phi_i = hho.edge_basis(iF, i);
						b(offset_F + i) = hho.integrate_over_edge(iF, [&](double x, double y){
								return nTF.dot(kappa(x,y) * grad_exact_solution(x,y)) * phi_i(x,y);
							});
					}
				}
			}
		}
	}
	
  return b;
}

double DiffusionEquation::H1Norm(HybridCore& hho, const Eigen::VectorXd XTF) {
  const auto mesh = hho.get_mesh_ptr();
	double value = 0.0;
  for (size_t iT = 0; iT < mesh-> n_cells(); iT++) {

		// Compute cell quadrature points and values of cell basis functions and gradients
		std::vector<HybridCore::qrule> quadT = hho.cell_qrule(iT, 2*hho.Ldeg());
		std::vector<Eigen::MatrixXd> dphiT_quadT = hho.grad_basis_quad(iT, quadT, hho.nlocal_cell_dofs());

		// CELL CONTRIBUTION
		//
		// Stiffness matrix
		//		(note: we compute also zeros first row and column...)
		Eigen::MatrixXd MG = hho.gram_matrix(dphiT_quadT, dphiT_quadT, hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs(), quadT, true);

		// Cell indices
		Eigen::VectorXd XTF_T = XTF.segment(iT*hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs());

		value += XTF_T.dot(MG * XTF_T);

		// Faces contributions
		
		for (size_t ilF=0; ilF < mesh->cell(iT)->n_edges() ; ilF++){
			size_t iF = mesh->cell(iT)->edge(ilF)->global_index();
  		auto hF = mesh->edge(iF)->measure();
			// Face quadrature points and values of cell and face basis functions at these points
			std::vector<HybridCore::qrule> quadF = hho.edge_qrule(iF, 2*hho.K()+2);
			std::vector<Eigen::VectorXd> phiT_quadF = hho.basis_quad('T', iT, quadF, hho.nlocal_cell_dofs());
			std::vector<Eigen::VectorXd> phiF_quadF = hho.basis_quad('F', iF, quadF, hho.nlocal_edge_dofs());

			// Face, cell and face-cell Gram matrix on F
			Eigen::MatrixXd MFF = hho.gram_matrix(phiF_quadF, phiF_quadF, hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs(), quadF, true);
			Eigen::MatrixXd MFT = hho.gram_matrix(phiF_quadF, phiT_quadF, hho.nlocal_edge_dofs(), hho.nlocal_cell_dofs(), quadF, false);
			Eigen::MatrixXd MTT_on_F = hho.gram_matrix(phiT_quadF, phiT_quadF, hho.nlocal_cell_dofs(), hho.nlocal_cell_dofs(), quadF, true);

			// Face indices
			Eigen::VectorXd XTF_F = XTF.segment(hho.ntotal_cell_dofs() + iF * hho.nlocal_edge_dofs(), hho.nlocal_edge_dofs());

			double valF = XTF_T.dot(MTT_on_F * XTF_T) - 2 * XTF_F.dot(MFT * XTF_T) + XTF_F.dot(MFF * XTF_F);

			value += valF / hF;
		}
  }
  return sqrt(value);

}


Eigen::VectorXd DiffusionEquation::VertexValues(HybridCore& hho, const Eigen::VectorXd XTF, const char from_dofs) {
  const auto mesh = hho.get_mesh_ptr();
	Eigen::VectorXd function_vertex = Eigen::VectorXd::Zero(mesh->n_vertices());

	if (from_dofs == 'T'){
		// compute from the cell polynomials, averaging all values coming from the cells around each vertex
		for (size_t iV = 0; iV < mesh->n_vertices(); iV++){
			auto xV = mesh->vertex(iV)->coords();
			auto cList = mesh->vertex(iV)->get_cells();
			auto weight = cList.size();
			for (size_t ilT = 0; ilT < weight; ilT++){
				auto temp = hho.evaluate_in_cell(XTF, cList[ilT]->global_index(), xV.x(), xV.y());
				function_vertex(iV) += temp;
			}
			function_vertex(iV) = function_vertex(iV)/(weight);
		}

	}else{
		// compute from the edge polynomials, averaging all values coming from the edges around each vertex
		for (size_t iV = 0; iV < mesh->n_vertices(); iV++){
			auto xV = mesh->vertex(iV)->coords();
			auto eList = mesh->vertex(iV)->get_edges();
			double weight = eList.size();
			for (size_t ilF = 0; ilF < weight; ilF++){
				auto temp = hho.evaluate_in_edge(XTF, eList[ilF]->global_index(), xV.x(), xV.y());
				function_vertex(iV) += temp;
			}
			function_vertex(iV) = function_vertex(iV)/(weight);
		}
	}

	return function_vertex;
}

double DiffusionEquation::get_assembly_time() const {
  return assembly_time;
}

double DiffusionEquation::get_solving_time() const {
  return solving_time;
}

double DiffusionEquation::get_itime(size_t idx) const {
  return itime[idx];
}

double DiffusionEquation::get_solving_error() const {
  return solving_error;
}

//@}
} // end of namespace HArDCore2D

#endif //_DIFFUSIONEQUATION_HPP
