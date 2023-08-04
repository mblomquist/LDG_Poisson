#include <iostream>
#include <functional>

#include "uvector.hpp"
#include "grid.hpp"
#include "Poisson.hpp"
#include "multigrid.hpp"

#define PI 3.14159265358979323846

using algoim::uvector;
using algoim::MultiLoop;

int main() {

    std::cout << "Hello, LDG Poisson Solver!" << std::endl;

    // Specify template parameters (order, dimensions)
    constexpr int P = 4;
    constexpr int N = 1;

//    std::cout << "\n--- Create a grid --- \n" << std::endl;
//    algoim::uvector<int, N> elements = 2;
//    algoim::uvector<double, N> domain_min = 0.;
//    algoim::uvector<double, N> domain_max = 1.;
//
//    PoissonSolver<P,N> solver;
//
//    solver.set_domain(domain_min, domain_max);
//    solver.set_elements_per_dim(elements);
//
//    std::cout << "\n--- Construct the Operators --- \n" << std::endl;
//    solver.construct_gradient_operator();
//
//    std::cout << "\n--- Project the rhs and true solution --- \n" << std::endl;
//    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return -4.*PI*PI*(sin(2.*PI*x(0)) + sin(2.*PI*x(1)));};
//    solver.project_rhs(cf_rhs);
//
//    std::function<double(uvector<double, N> x)> cf_sol = [](uvector<double, N> x) { return sin(2.*PI*x(0))+sin(2.*PI*x(1));};
//    solver.project_sol(cf_sol);
//
//    std::cout << "\n--- Printing Operators ---" << std::endl;
//    char output_operators[100];
//    sprintf(output_operators, "../out/operators.csv");
//    solver.print_operators_to_file(output_operators);
//
//    // print gradient operator in 3d
////    char output_G0[100];
////    char output_G1[100];
////    char output_G2[100];
////    sprintf(output_G0, "../out/G0.csv");
////    sprintf(output_G1, "../out/G1.csv");
////    sprintf(output_G2, "../out/G2.csv");
////    solver.print_3d_gradient_to_file(output_G0, output_G1, output_G2);
//
//    std::cout << "\n--- Printing Vectors ---" << std::endl;
//    char output_vectors[100];
//    sprintf(output_vectors, "../out/vectors.csv");
//    solver.print_vectors_to_file(output_vectors);

    test_transform<P,N>();

    return 0;
}
