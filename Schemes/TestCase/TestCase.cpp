// Class to Provides various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCase.hpp"
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <Eigen/Dense>


// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class
TestCase::TestCase(const std::vector<int> iTC)
	: iTC(iTC),
		_deg_diff(0) {
		validate();
		if (iTC[1]==2){
			_deg_diff = 2;
		}
	}

/////////////////// SOLUTION ///////////////////////////////:

// Solution
double TestCase::sol(const double x, const double y){
	double u = 0;
	switch(iTC[0]){
		case 1: u = sin(pi*x) * sin(pi*y);		/// iTC[0]=1: u(x,y)=sin(pi x) * sin(pi y)
						break;

		case 2: u = cos(pi*x) * cos(pi*y);		/// iTC[0]=2: u(x,y)=cos(pi x) * cos(pi y)
						break;

		case 3: u = x;		/// iTC[0]=3: u(x,y)= x
						break;

		case 4: u = y;		/// iTC[0]=3: u(x,y)= y
						break;

		case 5: u = pow(x, 2) + pow(y, 2); /// iTC[0]=4: u(x,y)= x^2 + y^2
						break;

		default: break;
	}
	return u;
}

// Gradient of the solution
Eigen::VectorXd TestCase::grad_sol(const double x, const double y){
	Eigen::VectorXd G = Eigen::VectorXd::Zero(2);
	switch(iTC[0]){
		case 1: G(0) = pi * cos(pi*x) * sin(pi*y);
						G(1) = pi * sin(pi*x) * cos(pi*y);
						break;

		case 2: G(0) = -pi * sin(pi*x) * cos(pi*y);
						G(1) = -pi * cos(pi*x) * sin(pi*y) ;
						break;

		case 3: G(0) = 1;
						G(1) = 0;
						break;

		case 4: G(0) = 0;
						G(1) = 1;
						break;

		case 5: G(0) = 2*x;
						G(1) = 2*y;
						break;

		default: break;
	}
	return G;
}

// Hessian of the solution
Eigen::MatrixXd TestCase::hess_sol(const double x, const double y){
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2,2);
	switch(iTC[0]){

		case 1: H.row(0) << - pi*pi*sin(pi*x)*sin(pi*y), pi*pi*cos(pi*x)*cos(pi*y);
						H.row(1) <<  pi*pi*cos(pi*x)*cos(pi*y), -pi*pi*sin(pi*x)*sin(pi*y);
						break;

		case 2: H.row(0) << - pi*pi*cos(pi*x)*cos(pi*y), pi*pi*sin(pi*x)*sin(pi*y);
						H.row(1) <<  pi*pi*sin(pi*x)*sin(pi*y), -pi*pi*cos(pi*x)*cos(pi*y);
						break;

		case 3: break;
		case 4: break;
		case 5: H.row(0) << 2, 0;
						H.row(1) << 0, 2;
						break;

		default: break;
	}
	return H;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
Eigen::MatrixXd TestCase::diff(const double x, const double y){
	Eigen::MatrixXd K = Eigen::MatrixXd::Identity(2,2);
	switch(iTC[1]){
		case 1: break;		/// iTC[1]=1: Diff = Id
		case 2: K.row(0) << pow(y,2)+1, -x*y;				/// iTC[1]=2: Diff = [y^2+1  -xy; -xy  x^2+1]
						K.row(1) << -x*y , pow(x,2)+1;
						break;
		case 3: if (x>0.5){				/// iTC[1]=3: Diff=Id if x<=1/2, Diff=[2 0; 0 1] if x>1/2. Only valid with iTC[0]=1;
							K.row(0) << 2, 0;
						}else{
							K.row(0) << 1, 0;
						}
						K.row(1) << 0, 1;
						break;
		default: break;
	}
	return K;
}

// Divergence by row of the diffusion matrix
Eigen::VectorXd TestCase::div_diff(const double x, const double y){
	Eigen::VectorXd divK = Eigen::VectorXd::Zero(2);
	switch(iTC[1]){
		case 1: break;
		case 2: divK(0) = -x;
						divK(1) = -y;
						break;
		case 3: break;
		default: break;
	}
	return divK;
}

///////////////////////////// SOURCE TERM ///////////////////////////

// Source term
double TestCase::source(const double x, const double y){

	Eigen::MatrixXd AHu = diff(x,y) * hess_sol(x,y);

	return -AHu.trace() - div_diff(x,y).dot(grad_sol(x,y));
}



///////////////////////////// VALIDATION ////////////////////////////

void TestCase::validate(){
	
	if (iTC[0]>5 || iTC[1]>3 || (iTC[1]==3 && iTC[0] !=1)){
		std::cout << "Incorrect choice of test cases: iTC= " << iTC[0] << ", " << iTC[1] << "\n";
		exit(EXIT_FAILURE);
	}	

}

