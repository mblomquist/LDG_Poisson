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
    algoim::uvector<int, N> elements;
    algoim::uvector<double, N> domain_min;
    algoim::uvector<double, N> domain_max;

    elements = 3;
//    elements(0) = 5;
//    elements(1) = 2;

    domain_min = -1.;
    domain_max = 1.;

    uniformGrid<N> myGrid(elements, domain_min, domain_max);
    std::cout << myGrid.get_total_faces() << std::endl;

//    // Initialize Poisson Solver
//    PoissonSolver<P,N> mySolver;
//
//    std::function<double(uvector<double, N> x)> cf_rhs = [](uvector<double, N> x) { return x(0)*x(1)*x(2);};
//    mySolver.project_rhs(cf_rhs, myGrid);
//
//    uvector<int, N> eval_grid = 11;
//    char output_file[100];
//    sprintf(output_file, "../out/projected_rhs.vtk");
//
//    mySolver.print_rhs_on_uniform_grid(eval_grid, myGrid, output_file);

    return 0;
}
