// Implementation of the HHO scheme in 2D for the diffusion equation
//
//   { -div(K \grad(u)) = f,       inside Omega
//   { K \grad(u) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
// 
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include <iostream>

#include <boost/program_options.hpp>
#include <boost/timer.hpp>

#include "mesh2d.hpp"
#include "import_mesh2d.hpp"
#include "mesh2d_builder.hpp"

#include "DiffusionEquation2d.hpp"
#include "hybridcore2d.hpp"
#include "TestCase/TestCase2d.hpp"
#include "vtu_writer2d.hpp"
#include "quad2d.hpp"

using namespace HArDCore2D;

// Mesh filenames
const std::string mesh_dir = "../../typ2_meshes/";
std::string default_mesh = mesh_dir + "cart5x5.typ2";
const std::string default_solver_type = "bicgstab";

constexpr int d = 2;

int main(int argc, const char* argv[]) {

  // --------------------------------------------------------------------------
  //                          Program options
  // --------------------------------------------------------------------------

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "Display this help message")
      ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh")
      ("edgedegree,k", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree on the edges")
      ("celldegree,l", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree in the cells")
      ("BC,b", boost::program_options::value<size_t>()->default_value(0), "Set the boundary conditions (0=Dirichlet, 1=Neumann)")
      ("testcase,c", boost::program_options::value<std::vector<int>>()->multitoken(), "Set the test case (as '-c n m'; n=exact sol, m=diffusion)")
      ("plot,p", boost::program_options::value<std::string>()->default_value("solution"), "Save plot of the solution to the given filename")
      ("solver_type", boost::program_options::value<std::string>()->default_value("bicgstab"), "Defines the linear solver for the system");


  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh
	std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);

  // Select the degree
  size_t K = vm["edgedegree"].as<size_t>();
  int L = vm["celldegree"].as<size_t>();

	// Select solver type
	std::string solver_type = (vm.count("solver_type") ? vm["solver_type"].as<std::string>() : default_solver_type);

	// Checks
	if ( (abs(int(K)-L) > 1) || (K<0) || (L<-1) ){
		std::cout << "Degrees k and l not in acceptable range (k positive and l=k-1, k, or k+1): k=" << K << ", l=" << L << "\n\n";
std::cout << std::abs(int(K)-L);
		exit(1);
	}

	// --------------------------------------------------------------------------
  //                        Create the HHO data structure
  // --------------------------------------------------------------------------

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
	if (mesh_ptr == NULL) {
		printf(
		  "Mesh cannot be created!\n Check the input file contains \n "
		  "Vertices "
		  "and cells with the correct tags");
		return 0;
	} 
 
  // Create the HHO structure
  HybridCore hho(mesh_ptr, K, L);

  // --------------------------------------------------------------------------
  //                        Create the model equation
  // --------------------------------------------------------------------------

  // Select the boundary conditions
	size_t BC = vm["BC"].as<size_t>();
	
	const std::vector<int> default_id_tcase = std::vector<int>{1,1};
	std::vector<int> id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_id_tcase);
	TestCase tcase(id_tcase);

  // Diffusion tensor
  DiffusionEquation::tensor_function_type kappa = [&](double x, double y) {
			return tcase.diff(x,y);
  };

  // Source term
  DiffusionEquation::scalar_function_type source = [&](double x, double y) {
			return tcase.source(x,y);
  };

  // Exact solution and gradient
  DiffusionEquation::scalar_function_type exact_solution = [&](double x, double y) {
			return tcase.sol(x,y);
  };
  DiffusionEquation::vector_function_type grad_exact_solution = [&](double x, double y) {
			return tcase.grad_sol(x,y);
  };

  // Create the model equation
  DiffusionEquation model(kappa, source, BC, exact_solution, grad_exact_solution, solver_type);


  // --------------------------------------------------------------------------
	// 											Recalling the mesh and parameters
  // --------------------------------------------------------------------------
	if (BC==0){
		std::cout << "\nBoundary conditions: Dirichlet\n";
	} else if (BC==1){
		std::cout << "\nBoundary conditions: Neumann\n";
	}
	std::cout << "Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "\n";
	size_t found = mesh_file.find_last_of("/\\");
	std::string mesh_name = mesh_file.substr(found+1);
	std::cout << "Mesh = " << mesh_name << " (nb cells= " << mesh_ptr->n_cells() << ", nb edges= " << mesh_ptr->n_edges() << ")\n";
	size_t nbedgedofs = mesh_ptr->n_edges()*hho.nlocal_edge_dofs();
	std::cout << "Degrees: edge = " << K << "; cell = " << L << " | Nb edge DOFs = " << nbedgedofs << "\n\n";

  // --------------------------------------------------------------------------
  //                        Solve the model problem
  // --------------------------------------------------------------------------

  std::cout << std::string(80, '-') << "\nAssembling and solving the problem..." << std::endl;
  boost::timer timer;
  Eigen::VectorXd u = model.solve(hho);
  std::cout<<"Model is solved"<<std::endl;
  
	std::cout << "Solver= " << solver_type << ", solved problem in " << timer.elapsed() << " seconds\n";
  std::cout << "\tAssembly time = " << model.get_assembly_time() << '\n';
  std::cout << "\tSolving time = " << model.get_solving_time() << '\n';
  std::cout << "\tResidual of the linear system = " << model.get_solving_error() << '\n';

	printf("interm times = %f | %f | %f | %f | %f\n", 
			model.get_itime(0),  model.get_itime(1),  model.get_itime(2),  model.get_itime(3), model.get_itime(4));

  // --------------------------------------------------------------------------
  //                        Compute the error
  // --------------------------------------------------------------------------

  std::cout << "Computing error..." << std::endl;

  // Interpolate the exact solution
  //timer.restart();
  Eigen::VectorXd UTF = hho.interpolate(exact_solution, 2*hho.K()+5);
  std::cout << "Interpolant calculated"<<std::endl;
	//std::cout << "Interpolant of exact solution computed in " << timer.elapsed() << " seconds.\n";

  // Compute the L2 error
  //timer.restart();
  double L2error = hho.L2norm(u - UTF)/hho.L2norm(UTF);
  //std::cout << "L2 Error = " << L2error << ", computed in " << timer.elapsed() << " seconds.\n";
	std::cout << "L2 Error = " << L2error << "\n";

  // Compute the L2 and H1 norms error
  //timer.restart();
  double H1error = model.H1Norm(hho,u - UTF)/model.H1Norm(hho,UTF);
  //std::cout << "H1 error = " << H1error << ", computed in " << timer.elapsed() << " seconds.\n";
	std::cout << "H1 Error = " << H1error << "\n";

  // Compute the Linf error of the edge coefficients
  //timer.restart();
//  double erroredges = hho.Linf_edge(u - UTF)/hho.Linf_edge(UTF);
//	std::cout << "Linf edges = " << erroredges << "\n";

  // Compute the integral
//  double average = hho.integral(u);
//	double average_exact_sol = hho.integrate_over_domain([&](auto x, auto y) {
//		return exact_solution(x,y);
//		});
//  std::cout << "Integral of solution = " << average << " (should be " << average_exact_sol << ")\n"; 

	// --------------------------------------------------------------------------
  //                     Creates a .vtu file of the solution
  // --------------------------------------------------------------------------

  if (vm.count("plot")) {
		std::string filename = vm["plot"].as<std::string>() + std::string(".vtu");
		VtuWriter2D plotdata(mesh_ptr);
		
		// Approximate solution, plots using cell or edge polynomials
		Eigen::VectorXd approx_sol_vertex_byT = model.VertexValues(hho, u, 'T');
		plotdata.write_to_vtu("T-" + filename,approx_sol_vertex_byT,1);

		Eigen::VectorXd approx_sol_vertex_byE = model.VertexValues(hho, u, 'E');
		plotdata.write_to_vtu("E-" + filename,approx_sol_vertex_byE,1);

		// Exact solution
		Eigen::VectorXd exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
		for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
			auto v= mesh_ptr->vertex(iV)->coords();
			exact_sol_vertex(iV) = exact_solution(v(0),v(1));
		}
		plotdata.write_to_vtu(std::string("exact-")+filename,exact_sol_vertex,1);

		std::cout << "\nApproximate solution in 'T-/E-" << filename << "', exact solution in '" << std::string("exact-")+filename << "'\n";

}

		// --------------------------------------------------------------------------
	  //                     Creates .txt file with data and results
	  // --------------------------------------------------------------------------

    std::ofstream out("results.txt");
    out << "BC: " << BC << '\n';
    out << "Solution: " << id_tcase[0] << '\n';
    out << "Diffusion: " << id_tcase[1] << '\n';
		out << "Mesh: " << mesh_name << "\n";
		out << "EdgeDegree: " << K << "\n";
		out << "CellDegree: " << L << "\n";
		out << "AssemblyTime: " << model.get_assembly_time() << "\n";
		out << "Solving time: " << model.get_solving_time() << "\n";
		out << "L2error: " << L2error << "\n";
		out << "H1error: " << H1error << "\n";
		out << "MeshSize: " << mesh_ptr->h_max() << "\n";
		out << "NbCells: " << mesh_ptr->n_cells() << "\n";
		out << "NbEdges: " << mesh_ptr->n_edges() << "\n";
		out << "NbEdgeDOFs: " << nbedgedofs << "\n";
		out << std::flush;
    out.close();


}
