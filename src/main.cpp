#include <iostream>
#include <functional>

#include "uvector.hpp"
#include "grid.hpp"
#include "Poisson.hpp"
#include "scratch.hpp"

#define PI 3.14159265358979323846

using algoim::uvector;
using algoim::MultiLoop;

int main() {

    std::cout << "Hello, LDG Poisson Solver!" << std::endl;

    // Specify template parameters (order, dimensions)
    constexpr int P = 3;
    constexpr int N = 2;

    std::cout << "\n--- Create a grid --- \n" << std::endl;
    algoim::uvector<int, N> elements = 2;
    algoim::uvector<double, N> domain_min = -1.;
    algoim::uvector<double, N> domain_max =  1.;

    PoissonSolver<P,N> solver;

    solver.set_domain(domain_min, domain_max);
    solver.set_elements_per_dim(elements);
    solver.construct_gradient_operator();
    solver.construct_penalty_operator();

    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return -PI*PI*(sin(PI*x(0))+sin(PI*x(1)));};
//    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return 1.;};
    solver.project_rhs(cf_rhs);

    std::function<double(uvector<double, N> x)> cf_sol = [](uvector<double, N> x) { return sin(PI*x(0))+sin(PI*x(1));};
    solver.project_sol(cf_sol);

    char output_operators[100];
    sprintf(output_operators, "../out/operators.csv");
    solver.print_operators_to_file(output_operators);

    char output_vectors[100];
    sprintf(output_vectors, "../out/vectors.csv");
    solver.print_vectors_to_file(output_vectors);

    return 0;
}
