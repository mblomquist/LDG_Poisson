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
    constexpr int levs = 4;
    algoim::uvector<int, N> elements = ipow(2,levs-1);
    algoim::uvector<double, N> domain_min = 0.;
    algoim::uvector<double, N> domain_max = 1.;
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
//    char output_G0[100];
//    char output_G1[100];
//    char output_G2[100];
//    sprintf(output_G0, "../out/G0.csv");
//    sprintf(output_G1, "../out/G1.csv");
//    sprintf(output_G2, "../out/G2.csv");
//    solver.print_3d_gradient_to_file(output_G0, output_G1, output_G2);
//
//    std::cout << "\n--- Printing Vectors ---" << std::endl;
//    char output_vectors[100];
//    sprintf(output_vectors, "../out/vectors.csv");
//    solver.print_vectors_to_file(output_vectors);

//    test_transform<P,N>();

//    constexpr int size = 3;
//    smatrix<double, size> A;
//    uvector<double, size> x, b;
//
//    A(0,0) = 2.; A(0,1) = -1.; A(0,2) = 0.;
//    A(1,0) = -1.; A(1,1) = 2.; A(1,2) = -1.;
//    A(2,0) = 0.; A(2,1) = -1.; A(2,2) = 2.;
//
//    b(0) = 0.9649; b(1) = 0.1576; b(2) = 0.9706;
//
//    x = 0.;
//    point_Gauss_Seidel<size>(A, x, b);
//
//    x = 0.;
//    double w = 1.1;
//    SOR<size>(A, x, b, w);

    uniformGrid<N> fineGrid(elements, domain_min, domain_max);

    BlockSparseMatrix<smatrix<double, ipow(P,N)>> T[N];

    MultiGrid<P,N> solver(fineGrid, T);

//    build_interpolation_operator<P,N>(fineGrid, I_cf);
//
//    for (int i = 0; i < fineGrid.get_total_elements(); ++i)
//    {
//        for (auto j : I_cf.row[i])
//        {
//            std::cout << "(" << i << ", " << j << ")" << std::endl;
//            I_cf(i,j).print();
//            std::cout << std::endl;
//        }
//    }
//
//    smatrix<double, 3> A, Ap;
//    A(0,0) = 1.; A(0, 1) = 2.; A(0, 2) = 3.;
//    A(1,0) = 2.; A(1, 1) = 4.; A(1, 2) = 6.;
//    A(2,0) = 3.; A(2, 1) = 6.; A(2, 2) = 9.;
//
//    Ap = pseudo_inverse_with_Eigen<3>(A);
//    Ap.print();

    return 0;
}
