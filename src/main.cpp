#include <iostream>
#include <functional>

#include "../../eigen/Eigen/Dense"
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
    constexpr int P = 2;
    constexpr int N = 2;

    std::cout << "\n--- Create a grid --- \n" << std::endl;
    constexpr int levs = 3;
    algoim::uvector<int, N> elements = ipow(2,levs-1);
    algoim::uvector<double, N> domain_min = 0.;
    algoim::uvector<double, N> domain_max = 1.;

    PoissonSolver<P,N> solver;

    solver.set_domain(domain_min, domain_max);
    solver.set_elements_per_dim(elements);

    std::cout << "\n--- Project the rhs and true solution --- \n" << std::endl;
    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return -4.*PI*PI*(sin(2.*PI*x(0)) + sin(2.*PI*x(1)));};
//    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return PI*PI*sin(PI*x(0));};
    solver.project_rhs(cf_rhs);

    std::function<double(uvector<double, N> x)> cf_sol = [](uvector<double, N> x) { return sin(2.*PI*x(0))+sin(2.*PI*x(1));};
//    std::function<double(uvector<double, N> x)> cf_sol = [](uvector<double, N> x) { return 0.;};
    solver.project_tsol(cf_sol);

    std::cout << "\n--- Solve --- \n" << std::endl;
    solver.solve_with_Multigrid();
    solver.compute_l2_error();

//    char op_file[100];
//    sprintf(op_file, "../out/op_file.csv");
//    solver.print_operators_to_file(op_file);
//
//    char v_file[100];
//    sprintf(v_file, "../out/v_file.csv");
//    solver.print_vectors_to_file(v_file);

    return 0;
}
