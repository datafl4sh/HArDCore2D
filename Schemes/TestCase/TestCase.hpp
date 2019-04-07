// Class to provide various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef _TEST_CASE_HPP
#define _TEST_CASE_HPP

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

/*!
* @defgroup TestCases
*	@brief Defines test cases (exact solutions, source terms etc.)
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

// @addtogroup TestCases
//@{

/// The TestCase class provides definition of test cases
class TestCase {

public:
  /// Initialise data
  TestCase(
    const std::vector<int> iTC  ///< The vector id of the test case: (id of solution, id of diffusion)
  );

	/// Returns the exact solution at the points x, y
	double sol(
		const double x,
		const double y
	);

	/// Returns the gradient of the exact solution at the points x, y
	Eigen::VectorXd grad_sol(
		const double x,
		const double y
	);

	/// Returns the Hessian of the exact solution at the points x, y
	Eigen::MatrixXd hess_sol(
		const double x,
		const double y
	);


	/// Returns the diffusion matrix at the points x, y
	Eigen::MatrixXd diff(
		const double x,
		const double y
	);

	/// Returns the divergence by row of the diffusion matrix at the points x, y
	Eigen::VectorXd div_diff(
		const double x,
		const double y
	);

	/// Returns the source term at the points x, y
	double source(
		const double x,
		const double y
	);

	/// Check if the provided test cases are valid (within range, and combination of solution/diffusion valid)
	void validate();

	/// Returns the degree of the diffusion tensor (useful to set up quadrature rules of proper degree)
	inline size_t get_deg_diff();

private:
  // Parameters: id of test case, pi
  const std::vector<int> iTC;
	const double pi = acos(-1);

	size_t _deg_diff;

};

inline size_t TestCase::get_deg_diff() { return _deg_diff; };

//@}

#endif //_TEST_CASE_HPP
