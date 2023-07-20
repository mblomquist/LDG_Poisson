#include <iostream>
#include <functional>

#include "uvector.hpp"
#include "grid.hpp"
#include "Poisson.hpp"

#define PI 3.14159265358979323846

using algoim::uvector;
using algoim::MultiLoop;

int main() {

    std::cout << "Hello, LDG Poisson Solver!" << std::endl;

    // Specify template parameters (order, dimensions)
    constexpr int P = 3;
    constexpr int N = 3;

    std::cout << "\n--- Create a grid --- \n" << std::endl;
    algoim::uvector<int, N> elements = 3;
    algoim::uvector<double, N> domain_min = -1.;
    algoim::uvector<double, N> domain_max =  1.;

    PoissonSolver<P,N> solver;

    solver.set_domain(domain_min, domain_max);
    solver.set_elements_per_dim(elements);

    solver.get_domain();

    char output_grid[100];
    sprintf(output_grid, "../out/grid.vtk");
    solver.print_grid(output_grid);


    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return x(0)*x(1)*x(2);};
    solver.project_rhs(cf_rhs);

    uvector<int, N> eval_grid = 11;

    char output_file[100];
    sprintf(output_file, "../out/rhs.vtk");
    solver.print_rhs_on_uniform_grid(eval_grid,  output_file);

    return 0;
}
