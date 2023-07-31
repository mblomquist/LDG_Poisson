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
    constexpr int P = 4;
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

    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return x(0)*x(1);};
    solver.project_rhs(cf_rhs);

    char output_rhs[100];
    sprintf(output_rhs, "../out/rhs.csv");
    solver.print_rhs_to_file(output_rhs);

    char output_grad[100];
    sprintf(output_grad, "../out/gradient.csv");
    solver.print_gradient_operator_to_file(output_grad);

    char output_penalty[100];
    sprintf(output_penalty, "../out/penalty.csv");
    solver.print_penalty_operator_to_file(output_penalty);

    return 0;
}
