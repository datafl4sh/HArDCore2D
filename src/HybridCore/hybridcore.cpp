// Core data structures and methods required to implement hybrid schemes in 2D (polynomial unknowns
// in the cells and on the edges, such as Hybrid High-order (HHO) schemes).
//
// Provides:
//  - Hybrid polynomial basis functions (on the cells and faces of the mesh)
//  - Generic routines to create quadrature nodes over cells and faces of the mesh
//  - Interpolation of general functions onto the HHO space
//  - Methods for integrating, evaluating, and computing norms of HHO solutions
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include <hybridcore.hpp>
#include <quad2d.hpp>
#include <quad1d.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <utils.hpp>
using namespace HArDCore2D;

// Creation class

HybridCore::HybridCore(const Mesh* mesh_ptr, const size_t K, const size_t L)
  : _mesh_ptr(mesh_ptr),
		_K(K),
		_L(L),
		_Ldeg(std::max(_L,0)),
		_nlocal_cell_dofs(dim_Pcell(_Ldeg)),
		_nlocal_edge_dofs(dim_Pedge(K)),
		_nhighorder_dofs(dim_Pcell(K+1)),
		_ngradient_dofs(_nhighorder_dofs - 1),
		_ntotal_cell_dofs(_nlocal_cell_dofs * mesh_ptr->n_cells()),
		_ntotal_edge_dofs(_nlocal_edge_dofs * mesh_ptr->n_edges()),
		_ninternal_edge_dofs(_nlocal_edge_dofs * mesh_ptr->n_i_edges()),
		_nboundary_edge_dofs(_nlocal_edge_dofs * mesh_ptr->n_b_edges()),
		_ntotal_dofs(_ntotal_cell_dofs+_ntotal_edge_dofs),
		_cell_bases(0),
		_cell_gradients(0),
		_edge_bases(0)	  {
	    std::cout << "Construct HybridCore" << std::endl;
			// Initialise basis functions on cells
			for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++) {
				auto cell_basis = create_cell_basis(iT);
				_cell_bases.push_back(cell_basis.first);
				_cell_gradients.push_back(cell_basis.second);
			}
			// Initialse basis fucntion on edges
			for (size_t iF = 0; iF < mesh_ptr->n_edges(); iF++){
				_edge_bases.push_back(create_edge_basis(iF));
			}
}

// -------------------------------------------------
// ------- Cell and edge basis functions
// -------------------------------------------------

size_t HybridCore::dim_Pcell(const size_t m) const {
	return (m + 1) * (m + 2) / 2;
}

size_t HybridCore::dim_Pedge(const size_t m) const {
	return m+1;
}

const HybridCore::cell_basis_type& HybridCore::cell_basis(size_t iT, size_t i) const {
  assert(iT < _mesh_ptr->n_cells());
  assert(i < _cell_bases[iT].size());
  return _cell_bases[iT][i];
}

const HybridCore::edge_basis_type& HybridCore::edge_basis(size_t iF, size_t i) const {
  assert(iF < _mesh_ptr->n_edges());
  assert(i < _edge_bases[iF].size());
  return _edge_bases[iF][i];
}

const HybridCore::cell_gradient_type& HybridCore::cell_gradient(size_t iT,
                                                                size_t i) const {
  assert(iT < _mesh_ptr->n_cells());
  assert(i < _cell_gradients[iT].size());
  return _cell_gradients[iT][i];
}
// const Mesh& HybridCore::mesh() const { return _mesh; }
std::pair<std::vector<HybridCore::cell_basis_type>,
			std::vector<HybridCore::cell_gradient_type> >
	HybridCore::create_cell_basis(size_t iT) const {
		std::vector<cell_basis_type> cell_basis;
		std::vector<cell_gradient_type> cell_gradient;
		auto xT = _mesh_ptr->cell(iT)->center_mass();
		auto hT = _mesh_ptr->cell(iT)->diam();

		// Create basis functions up to degree K+1 since we need a high-order
		// basis for the reconstruction operators
		for (size_t degree = 0; degree <= _K+1; degree++){
			for (size_t i = 0; i <= degree; i++) {
				size_t j = degree - i;
				cell_basis_type phi = [ xT, i, j, hT ](double x, double y)->double {
				      return std::pow( (x - xT.x())/hT , i) * std::pow( (y - xT.y())/hT , j);
			  		};

			  cell_gradient_type dphi = [ xT, hT, i, j ](double x, double y)->Eigen::Vector2d {
			      Vector2d gradient;
			      gradient(0) =
									i == 0 ? 0.0 : i * std::pow( (x - xT.x())/hT , i - 1) *
						                             std::pow( (y - xT.y())/hT , j) / hT;
			      gradient(1) =
				          j == 0 ? 0.0 : std::pow( (x - xT.x())/hT , i) * j *
						                             std::pow( (y - xT.y())/hT , j - 1) / hT;
						return gradient;
					};

					cell_basis.push_back(std::move(phi));
					cell_gradient.push_back(std::move(dphi));		 
			}
		}
		return std::make_pair(std::move(cell_basis), std::move(cell_gradient));
}

std::vector<HybridCore::edge_basis_type> HybridCore::create_edge_basis(size_t iF) const{
	std::vector<edge_basis_type> edge_basis;
	
	auto xF = _mesh_ptr->edge(iF)->center_mass();
	auto edge_tang = _mesh_ptr->edge(iF)->tangent();
	auto hF = _mesh_ptr->edge(iF)->measure();
	
	edge_tang = edge_tang.normalized();
	for (size_t degree = 0; degree <= _K; degree++){
		edge_basis_type phi = [xF, edge_tang, hF, degree](double x, double y)
					->double{
						Eigen::Vector2d edge_vec = Eigen::Vector2d(x,y) - xF;
						double s = edge_vec.dot(edge_tang);
						return std::pow( s/hF , degree);
					};
					edge_basis.push_back(std::move(phi));
	}
	return std::move(edge_basis);
}

// -------------------------------------------------
// --------- Create quadrature rules
// -------------------------------------------------


std::vector<HybridCore::qrule> HybridCore::cell_qrule(const size_t iT,
                                                      const size_t doe) const {
  assert(iT < _mesh_ptr->n_cells());

  QuadRuleTriangle quadCell(std::max(int(doe),1) + _offset_doe, true);
  std::vector<qrule> quad;
  const Cell* cell = _mesh_ptr->cell(iT);
  size_t nedges = cell->n_edges();

  if (nedges == 3) {
      // Triangle
      auto x0 = cell->vertex(0)->coords();
      auto x1 = cell->vertex(1)->coords();
      auto x2 = cell->vertex(2)->coords();

      double xT[] = {x0.x(), x1.x(), x2.x()};
      double yT[] = {x0.y(), x1.y(), x2.y()};

      quadCell.setup(xT, yT);
      for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
          quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.wq(iqn));
      }

  } else if (nedges == 4) {  
			// rect split into two tgls
      auto x0 = cell->vertex(0)->coords();
      for (size_t isplit = 0; isplit < 2; isplit++) {
          auto x1 = cell->vertex(1 + isplit)->coords();
          auto x2 = cell->vertex(2 + isplit)->coords();
          double xTr[] = {x0.x(), x1.x(), x2.x()};
          double yTr[] = {x0.y(), x1.y(), x2.y()};
          quadCell.setup(xTr, yTr);

          for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
              quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.wq(iqn));
          }
      }


  } else {
      // other split at barycentre size_t smaller triangles
      auto xF = cell->center_mass();

      for (size_t ilF = 0; ilF < nedges; ilF++) {
          const Edge* e = cell->edge(ilF);
          auto x0 = e->vertex(0)->coords();
          auto x1 = e->vertex(1)->coords();

          double xTr[] = {xF.x(), x0.x(), x1.x()};
          double yTr[] = {xF.y(), x0.y(), x1.y()};
          quadCell.setup(xTr, yTr);

          for (size_t iqn = 0; iqn < quadCell.nq(); iqn++) {
              quad.emplace_back(quadCell.xq(iqn), quadCell.yq(iqn), quadCell.wq(iqn));
          }
      }
  }

  return quad;
}

std::vector<HybridCore::qrule> HybridCore::edge_qrule(const size_t iE,
                                                      const size_t doe) const {
    assert(iE < _mesh_ptr->n_edges());
    QuadRuleEdge quadEdge(std::max(int(doe),1) + _offset_doe, true);
    std::vector<qrule> quad;
    Edge* iedge = _mesh_ptr->edge(iE);

    auto x0 = iedge->vertex(0)->coords();
    auto x1 = iedge->vertex(1)->coords();
    double xT[] = {x0.x(), x1.x()};
    double yT[] = {x0.y(), x1.y()};
	
    quadEdge.setup(xT, yT);
    for (size_t iqn = 0; iqn < quadEdge.nq(); iqn++) {
        quad.emplace_back(quadEdge.xq(iqn), quadEdge.yq(iqn), quadEdge.wq(iqn));
    }
    return quad;
}

// -------------------------------------------------------------------
// ------  Gram matrices for scalar and vector-valued functions
// -------------------------------------------------------------------

Eigen::MatrixXd HybridCore::gram_matrix(const std::vector<Eigen::ArrayXd>& f_quad, const std::vector<Eigen::ArrayXd>& g_quad, const size_t& nrows, const size_t& ncols, const std::vector<HybridCore::qrule>& quad, const bool& sym, std::vector<double> L2weight) const {

	Eigen::MatrixXd GGM = Eigen::MatrixXd::Zero(nrows, ncols);

	size_t nbq = quad.size();

	// Recast product of quadrature and L2weight into an Eigen::ArrayXd
	Eigen::ArrayXd quad_L2_weights = Eigen::ArrayXd::Zero(nbq);
	if (L2weight.size() == 0){
		for (size_t iqn = 0; iqn < nbq; iqn++){
			quad_L2_weights(iqn) = quad[iqn].w;
		}
	}else{
		for (size_t iqn = 0; iqn < nbq; iqn++){
			quad_L2_weights(iqn) = quad[iqn].w * L2weight[iqn];
		}
	}

	for (size_t i = 0; i < nrows; i++){
		size_t jcut = 0;
		if (sym) jcut = i;
		for (size_t j = 0; j < jcut; j++){
				GGM(i, j) = GGM(j, i);
		}
		for (size_t j = jcut; j < ncols; j++){
			// Integrate f_i * g_j
			// The products here are component-wise since the terms are Eigen::ArrayXd
			GGM(i, j) = (quad_L2_weights * f_quad[i] * g_quad[j]).sum();
		}
	}

	return GGM;
}

Eigen::MatrixXd HybridCore::gram_matrix(const std::vector<Eigen::ArrayXXd>& F_quad, const std::vector<Eigen::ArrayXXd>& G_quad, const size_t& nrows, const size_t& ncols, const std::vector<HybridCore::qrule>& quad, const bool& sym, std::vector<Eigen::Matrix2d> L2Weight) const {

	Eigen::MatrixXd GSM = Eigen::MatrixXd::Zero(nrows, ncols);

	if (L2Weight.size() == 0){
		L2Weight.resize(quad.size(), Eigen::Matrix2d::Identity());
	}

	size_t nbq = quad.size();
	for (size_t i = 0; i < nrows; i++){
		size_t jcut = 0;
		if (sym) jcut = i;
		for (size_t j = 0; j < jcut; j++){
				GSM(i, j) = GSM(j, i);
		}
		for (size_t j = jcut; j < ncols; j++){
			// Integrate f_i * g_j
			for (size_t iqn = 0; iqn < nbq; iqn++){
				GSM(i, j) += quad[iqn].w * (L2Weight[iqn] * F_quad[i].col(iqn).matrix()).dot(G_quad[j].col(iqn).matrix());
			}
		}
	}

	return GSM;
}

// ----------------------------------------------------------------
// ------- Basis functions and gradients on quadrature nodes
// ----------------------------------------------------------------

std::vector<Eigen::ArrayXd> HybridCore::basis_quad(const char B, const size_t iTF, const std::vector<qrule> quad, const size_t nbasis) const {
	size_t nbq = quad.size();
	std::vector<Eigen::ArrayXd> phi_quad(nbasis, Eigen::ArrayXd::Zero(nbq));	
	for (size_t i = 0; i < nbasis; i++){
		if (B == 'T'){
			const auto &phi_i = cell_basis(iTF, i);
			for (size_t iqn = 0; iqn < nbq; iqn++){
				phi_quad[i](iqn) = phi_i(quad[iqn].x, quad[iqn].y);
			}

		}else{
			const auto &phi_i = edge_basis(iTF, i);
			for (size_t iqn = 0; iqn < nbq; iqn++){
				phi_quad[i](iqn) = phi_i(quad[iqn].x, quad[iqn].y);
			}
		}
	}

	return phi_quad;
}

std::vector<Eigen::ArrayXXd> HybridCore::grad_basis_quad(const size_t iT, const std::vector<qrule> quad, const size_t nbasis) const {
	size_t nbq = quad.size();
	std::vector<Eigen::ArrayXXd> dphi_quad(nbasis, Eigen::ArrayXXd::Zero(_mesh_ptr->dim(), nbq));	

	// No need to consider i=0 since the first basis function is constant
	for (size_t i = 1; i < nbasis; i++){
		auto &dphi_i = cell_gradient(iT, i);

		for (size_t iqn=0; iqn<nbq; iqn++){
			dphi_quad[i].col(iqn) = dphi_i(quad[iqn].x, quad[iqn].y);
		}
	}		

	return dphi_quad;
}

// -----------------------------------------------------------------
// ------- Weights to eliminate cell unknows if L=-1
// -----------------------------------------------------------------


Eigen::VectorXd HybridCore::compute_weights(size_t iT) const {
	Cell* iCell = _mesh_ptr->cell(iT);
	size_t nlocal_edges = iCell->n_edges();
	Eigen::VectorXd barycoefT = Eigen::VectorXd::Zero(nlocal_edges);

	// Rule degree 0: all coefficients identical
//	barycoefT = (1/double(nlocal_faces)) * Eigen::VectorXd::Ones(nlocal_faces);

	// Rule degree 1: coefficient is |F|d_{TF}/(d|T|)
	for (size_t ilF=0; ilF < nlocal_edges; ilF++){
		Eigen::Vector2d normalTF = iCell->edge_normal(ilF);
		Eigen::Vector2d xFxT = iCell->edge(ilF)->center_mass() - iCell->center_mass();

		double dTF = xFxT.dot(normalTF); 

		barycoefT(ilF) = iCell->edge(ilF)->measure() * dTF / (_mesh_ptr->dim() * iCell->measure());

	}

	return barycoefT;
}

// ----------------------------------------------------------
// ---------- Norms of discrete unknowns
// ----------------------------------------------------------

Eigen::VectorXd HybridCore::restr(const Eigen::VectorXd &Xh, size_t iT) const {
	
	Cell* cell = _mesh_ptr->cell(iT);
	size_t nedgesT = cell->n_edges();

	Eigen::VectorXd XTF = Eigen::VectorXd::Zero(_nlocal_cell_dofs + nedgesT * _nlocal_edge_dofs);

	XTF.head(_nlocal_cell_dofs) = Xh.segment(iT * _nlocal_cell_dofs, _nlocal_cell_dofs);
	for (size_t ilF = 0; ilF < nedgesT; ilF++){
		size_t offset_F = _ntotal_cell_dofs + cell->edge(ilF)->global_index() * _nlocal_edge_dofs;
		XTF.segment(_nlocal_cell_dofs + ilF * _nlocal_edge_dofs, _nlocal_edge_dofs) =
				Xh.segment(offset_F, _nlocal_edge_dofs);
	}

	return XTF;
}

double HybridCore::L2norm(const Eigen::VectorXd &Xh) const {
  double value = 0.0;
  for (size_t iT = 0; iT < _mesh_ptr->n_cells(); iT++) {
		// L2 norm computed using the mass matrix
		// Compute cell quadrature nodes and values of cell basis functions at these nodes
		std::vector<HybridCore::qrule> quadT = cell_qrule(iT, 2*(_K+1));
		std::vector<Eigen::ArrayXd> phi_quadT = basis_quad('T', iT, quadT, _nlocal_cell_dofs);

		Eigen::MatrixXd MTT = gram_matrix(phi_quadT, phi_quadT, _nlocal_cell_dofs, _nlocal_cell_dofs, quadT, true);
		Eigen::VectorXd XT = Xh.segment(iT*_nlocal_cell_dofs,_nlocal_cell_dofs);

		value += XT.dot(MTT*XT);

  }
  return sqrt(value);
}

double HybridCore::H1norm(const Eigen::VectorXd &Xh) const {

	double value = 0.0;
  for (size_t iT = 0; iT < _mesh_ptr-> n_cells(); iT++) {
		Cell* cell = _mesh_ptr->cell(iT);
		size_t nedgesT = cell->n_edges();
		size_t nlocal_dofs = _nlocal_cell_dofs + nedgesT * _nlocal_edge_dofs;

		// Local matrix of the bilinear form corresponding to the cell contribution in the discrete H1 norm
		Eigen::MatrixXd H1aT = Eigen::MatrixXd::Zero(nlocal_dofs, nlocal_dofs);

		// Compute cell quadrature nodes and values of gradients of cell basis functions
		std::vector<HybridCore::qrule> quadT = cell_qrule(iT, 2*_Ldeg);
		std::vector<Eigen::ArrayXXd> dphiT_quadT = grad_basis_quad(iT, quadT, _nlocal_cell_dofs);

		// CELL CONTRIBUTION
		//
		// Stiffness matrix
		H1aT.topLeftCorner(_nlocal_cell_dofs, _nlocal_cell_dofs) = 
				gram_matrix(dphiT_quadT, dphiT_quadT, _nlocal_cell_dofs, _nlocal_cell_dofs, quadT, true);


		// FACES CONTRIBUTION		
		for (size_t ilF=0; ilF < nedgesT ; ilF++){
			size_t iF = cell->edge(ilF)->global_index();
			size_t offset_F = _nlocal_cell_dofs + ilF * _nlocal_edge_dofs;
  		auto hF = _mesh_ptr->edge(iF)->measure();
			// Face quadrature nodes and values of cell and face basis functions at these nodes
			std::vector<HybridCore::qrule> quadF = edge_qrule(iF, 2*_K+1);
			std::vector<Eigen::ArrayXd> phiT_quadF = basis_quad('T', iT, quadF, _nlocal_cell_dofs);
			std::vector<Eigen::ArrayXd> phiF_quadF = basis_quad('F', iF, quadF, _nlocal_edge_dofs);

			// Face, face-cell and face-face Gram matrix on F
			Eigen::MatrixXd MFF = gram_matrix(phiF_quadF, phiF_quadF, _nlocal_edge_dofs, _nlocal_edge_dofs, quadF, true);
			Eigen::MatrixXd MFT = gram_matrix(phiF_quadF, phiT_quadF, _nlocal_edge_dofs, _nlocal_cell_dofs, quadF, false);
			Eigen::MatrixXd MTT_on_F = gram_matrix(phiT_quadF, phiT_quadF, _nlocal_cell_dofs, _nlocal_cell_dofs, quadF, true);

			// Contribution of the face to the local bilinear form
			H1aT.block(offset_F, offset_F, _nlocal_edge_dofs, _nlocal_edge_dofs) += MFF / hF;
			H1aT.block(offset_F, 0, _nlocal_edge_dofs, _nlocal_cell_dofs) -= MFT / hF;
			H1aT.block(0, offset_F, _nlocal_cell_dofs, _nlocal_edge_dofs) -= MFT.transpose() / hF;
			H1aT.topLeftCorner(_nlocal_cell_dofs, _nlocal_cell_dofs) += MTT_on_F / hF;
		}

		Eigen::VectorXd XTF = restr(Xh, iT);
		value += XTF.transpose() * H1aT * XTF;
  }

  return sqrt(value);
}


double HybridCore::Linf_edge(const Eigen::VectorXd &Xh) const {
  double value = 0.0;

	Eigen::VectorXd XF = Xh.segment(_mesh_ptr->n_cells()*_nlocal_cell_dofs,_mesh_ptr->n_edges()*_nlocal_edge_dofs);

  for (size_t iF = _mesh_ptr->n_cells()*_nlocal_cell_dofs; iF < _mesh_ptr->n_edges()*_nlocal_edge_dofs; iF++) {
		value = std::max(value, std::abs(XF(iF)));

  }
  return value;
}




// -----------------------------------------------------------
// --------- Evaluate discrete functions in cell or edges
// -----------------------------------------------------------

double HybridCore::evaluate_in_cell(const Eigen::VectorXd XTF, size_t iT, double x, double y) const {
  double value = 0.0;
  for (size_t i = 0; i < _nlocal_cell_dofs; i++) {
    const auto &phi_i = cell_basis(iT, i);
    const size_t index = i + iT * _nlocal_cell_dofs;
    value += XTF(index) * phi_i(x,y);
  }
  return value;
}

double HybridCore::evaluate_in_edge(const Eigen::VectorXd XTF, size_t iF, double x, double y) const {
  double value = 0.0;
  for (size_t i = 0; i < _nlocal_edge_dofs; i++) {
    const auto &phi_i = edge_basis(iF, i);
    const size_t index = _ntotal_cell_dofs + iF * _nlocal_edge_dofs + i;
    value += XTF(index) * phi_i(x,y);
  }
  return value;
}


double HybridCore::integral(const Eigen::VectorXd &XTF) const {
  double value = 0.0;
  for (size_t iT = 0; iT < _mesh_ptr->n_cells(); iT++) {
    value += integrate_over_cell(iT, [&XTF, iT, this](double x, double y) {
      return evaluate_in_cell(XTF, iT, x, y);
    });
  }
  return value;
}


// ---------------------------------------------------
// ------- Functions that return class elements
// ---------------------------------------------------

const Mesh* HybridCore::get_mesh_ptr() {return _mesh_ptr;}
const size_t HybridCore::K() {return _K;}
const int HybridCore::L() {return _L;}
const int HybridCore::Ldeg() {return _Ldeg;}
