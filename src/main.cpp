#include <iostream>
#include <functional>

#include "uvector.hpp"
#include "multiloop.hpp"
#include "smatrix.hpp"
#include "legendre.hpp"
#include "grid.hpp"

#define PI 3.14159265358979323846

using algoim::uvector;
using algoim::MultiLoop;

int main() {

    std::cout << "Hello, World!" << std::endl;

    // Check templated legendre
    constexpr int P = 3;
    constexpr int N = 3;

    std::cout << "\n--- L2 Projection --- \n" << std::endl;
    std::function<double(uvector<double, N> x)> rhs = [](uvector<double, N> x) { return x(0)*x(1)*x(2);};

    std::cout << "\n--- Make a Grid --- \n" << std::endl;
    algoim::uvector<int, N> elements;
    algoim::uvector<double, N> domain_min;
    algoim::uvector<double, N> domain_max;

    elements(0) = 5;
    elements(1) = 2;
    elements(2) = 4;

    domain_min = -1.;
    domain_max = 1.;

    uniformGrid<N> myGrid(elements, domain_min, domain_max);
    char output_grid[100];
    sprintf(output_grid, "../out/grid.vtk");
    myGrid.print_grid_to_vtk(output_grid);

    std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> projected_func;
    project_func<P,N>(projected_func, rhs, myGrid);

    uvector<int, N> eval_grid = 11;
    char output_file[100];
    sprintf(output_file, "../out/projected_func.vtk");

    evaluate_basis_on_uniform_grid<P,N>(projected_func, myGrid, eval_grid, output_file, "point_data");

    return 0;
}
