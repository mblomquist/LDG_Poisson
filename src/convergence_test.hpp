#pragma once

/*
 * Order-of-accuracy benchmark for the LDG Poisson solver.
 *
 * Two families of manufactured solutions are used:
 *
 *   Periodic case:
 *       u(x) = prod_d sin(2*pi*x_d),   Laplacian(u) = -N*(2*pi)^2 * u,
 *   which is smooth and periodic on the unit cube and has zero mean (so it
 *   satisfies the compatibility condition of the singular periodic operator).
 *   This case exercises the original multigrid solver path unchanged.
 *
 *   Dirichlet / Neumann / Robin cases:
 *       u(x) = prod_d sin(pi*x_d),     Laplacian(u) = -N*pi^2 * u,
 *   which vanishes on the boundary of the unit cube (so homogeneous Dirichlet
 *   data is exact) and has a non-trivial normal derivative for the Neumann and
 *   Robin tests.  These cases use the direct sparse solver path.
 *
 * The solver convention is sol -> -u (compute_l2_error measures ||sol + t_sol||
 * with t_sol projected to +u).  Boundary data is therefore prescribed for the
 * field phi = sol = -u:
 *
 *     Dirichlet:  g = phi                       = -u
 *     Neumann:    g = d(phi)/dn                 = -n * du/dx_dim
 *     Robin:      g = alpha*phi + beta*d(phi)/dn = -(alpha*u + beta*n*du/dx_dim)
 *
 * The pure Neumann problem is singular; the solver appends a mean-zero
 * constraint and the error is measured after removing the constant mode from
 * both fields (subtract_mean = true).
 */

#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <functional>

#include "uvector.hpp"
#include "math_tools.hpp"
#include "Poisson.hpp"

namespace convergence_test
{
    constexpr double kPi = 3.14159265358979323846;

    // ---- Periodic manufactured solution: u(x) = prod_d sin(2*pi*x_d). --------
    template<int N>
    double manufactured_solution(const algoim::uvector<double, N> &x)
    {
        double u = 1.;
        for (int dim = 0; dim < N; ++dim) {
            u *= std::sin(2. * kPi * x(dim));
        }
        return u;
    }

    template<int N>
    double manufactured_rhs(const algoim::uvector<double, N> &x)
    {
        return -double(N) * (2. * kPi) * (2. * kPi) * manufactured_solution<N>(x);
    }

    // ---- Non-periodic manufactured solution: u(x) = prod_d sin(pi*x_d). ------
    template<int N>
    double bc_solution(const algoim::uvector<double, N> &x)
    {
        double u = 1.;
        for (int dim = 0; dim < N; ++dim) {
            u *= std::sin(kPi * x(dim));
        }
        return u;
    }

    // Laplacian(u) = -N*pi^2 * u  (the solver's right-hand side, sol -> -u).
    template<int N>
    double bc_rhs(const algoim::uvector<double, N> &x)
    {
        return -double(N) * kPi * kPi * bc_solution<N>(x);
    }

    // Partial derivative du/dx_dim of the non-periodic manufactured solution.
    template<int N>
    double bc_grad(int dim, const algoim::uvector<double, N> &x)
    {
        double g = kPi * std::cos(kPi * x(dim));
        for (int d = 0; d < N; ++d) {
            if (d != dim)
                g *= std::sin(kPi * x(d));
        }
        return g;
    }

    // Observed order of accuracy between two consecutive refinement levels.
    inline double observed_order(double h_coarse, double err_coarse,
                                 double h_fine, double err_fine)
    {
        return std::log(err_coarse / err_fine) / std::log(h_coarse / h_fine);
    }

    // ---- Periodic study (original multigrid solver path). --------------------
    template<int P, int N>
    void run_periodic_study(const std::vector<int> &resolutions,
                            std::vector<double> &hs,
                            std::vector<double> &errors)
    {
        hs.clear();
        errors.clear();

        for (int n_elem : resolutions) {
            std::cout << "\n--- periodic, " << N << "D, P = " << P
                      << ", elements/dim = " << n_elem << " ---" << std::endl;

            PoissonSolver<P, N> solver;

            algoim::uvector<int, N>    elements   = n_elem;
            algoim::uvector<double, N> domain_min = 0.;
            algoim::uvector<double, N> domain_max = 1.;

            solver.set_domain(domain_min, domain_max);
            solver.set_elements_per_dim(elements);

            std::function<double(const algoim::uvector<double, N> &)> cf_rhs =
                    [](const algoim::uvector<double, N> &x) { return manufactured_rhs<N>(x); };
            std::function<double(const algoim::uvector<double, N> &)> cf_sol =
                    [](const algoim::uvector<double, N> &x) { return manufactured_solution<N>(x); };

            solver.project_rhs(cf_rhs);
            solver.project_tsol(cf_sol);

            solver.solve_with_Multigrid();

            hs.push_back(1. / double(n_elem));
            errors.push_back(solver.compute_l2_error());
        }
    }

    // ---- Boundary-condition study (direct solver path). ----------------------
    // bc_type is one of BC_DIRICHLET, BC_NEUMANN, BC_ROBIN.
    template<int P, int N>
    void run_bc_study(int bc_type,
                      const std::vector<int> &resolutions,
                      std::vector<double> &hs,
                      std::vector<double> &errors,
                      double alpha = 1., double beta = 1.)
    {
        hs.clear();
        errors.clear();

        const bool subtract_mean = (bc_type == BC_NEUMANN);

        for (int n_elem : resolutions) {
            std::cout << "\n--- " << N << "D, P = " << P
                      << ", elements/dim = " << n_elem << " ---" << std::endl;

            PoissonSolver<P, N> solver;

            algoim::uvector<int, N>    elements   = n_elem;
            algoim::uvector<double, N> domain_min = 0.;
            algoim::uvector<double, N> domain_max = 1.;

            solver.set_domain(domain_min, domain_max);
            solver.set_elements_per_dim(elements);

            std::function<double(const algoim::uvector<double, N> &)> cf_rhs =
                    [](const algoim::uvector<double, N> &x) { return bc_rhs<N>(x); };
            std::function<double(const algoim::uvector<double, N> &)> cf_sol =
                    [](const algoim::uvector<double, N> &x) { return bc_solution<N>(x); };

            solver.project_rhs(cf_rhs);
            solver.project_tsol(cf_sol);

            // Prescribe the boundary data for the field phi = sol = -u.
            for (int dim = 0; dim < N; ++dim) {
                for (int side = 0; side < 2; ++side) {
                    double n = (side == 1) ? 1. : -1.;

                    if (bc_type == BC_DIRICHLET) {
                        std::function<double(const algoim::uvector<double, N> &)> g =
                                [](const algoim::uvector<double, N> &x) { return -bc_solution<N>(x); };
                        solver.set_bc(dim, side, BC_DIRICHLET, g);
                    }
                    else if (bc_type == BC_NEUMANN) {
                        std::function<double(const algoim::uvector<double, N> &)> g =
                                [dim, n](const algoim::uvector<double, N> &x) {
                                    return -n * bc_grad<N>(dim, x);
                                };
                        solver.set_bc(dim, side, BC_NEUMANN, g);
                    }
                    else { // BC_ROBIN
                        std::function<double(const algoim::uvector<double, N> &)> g =
                                [dim, n, alpha, beta](const algoim::uvector<double, N> &x) {
                                    return -(alpha * bc_solution<N>(x) + beta * n * bc_grad<N>(dim, x));
                                };
                        solver.set_bc(dim, side, BC_ROBIN, g, alpha, beta);
                    }
                }
            }

            solver.solve();

            hs.push_back(1. / double(n_elem));
            errors.push_back(solver.compute_l2_error(subtract_mean));
        }
    }

    // Print a table of (N, h, L2 error, observed order) for one case.
    inline void report_table(const std::string &label, int P,
                             const std::vector<int> &resolutions,
                             const std::vector<double> &hs,
                             const std::vector<double> &errors)
    {
        std::cout << "\n============================================================\n";
        std::cout << " " << label << " (P = " << P
                  << ", degree " << P - 1 << ", expected L2 order ~ " << P << ")\n";
        std::cout << "============================================================\n";
        std::cout << std::setw(8)  << "N"
                  << std::setw(16) << "h"
                  << std::setw(18) << "L2 error"
                  << std::setw(14) << "order" << "\n";
        std::cout << "------------------------------------------------------------\n";

        for (std::size_t i = 0; i < resolutions.size(); ++i) {
            std::cout << std::setw(8) << resolutions[i]
                      << std::setw(16) << std::scientific << std::setprecision(4) << hs[i]
                      << std::setw(18) << std::scientific << std::setprecision(6) << errors[i];

            if (i == 0) {
                std::cout << std::setw(14) << "-";
            } else {
                double order = observed_order(hs[i - 1], errors[i - 1], hs[i], errors[i]);
                std::cout << std::setw(14) << std::fixed << std::setprecision(3) << order;
            }
            std::cout << "\n";
        }
        std::cout << "============================================================\n";
    }

    // Run the full 2D convergence benchmark for all four boundary-condition
    // cases and report the tables.
    template<int P>
    void run_benchmark()
    {
        std::cout << "\n########## LDG Poisson convergence benchmark ##########\n";

        std::vector<int> res = {4, 8, 16, 32};

        std::vector<double> h_per, err_per;
        run_periodic_study<P, 2>(res, h_per, err_per);

        std::vector<double> h_dir, err_dir;
        run_bc_study<P, 2>(BC_DIRICHLET, res, h_dir, err_dir);

        std::vector<double> h_neu, err_neu;
        run_bc_study<P, 2>(BC_NEUMANN, res, h_neu, err_neu);

        std::vector<double> h_rob, err_rob;
        run_bc_study<P, 2>(BC_ROBIN, res, h_rob, err_rob, /*alpha=*/1., /*beta=*/1.);

        report_table("2D periodic",  P, res, h_per, err_per);
        report_table("2D Dirichlet", P, res, h_dir, err_dir);
        report_table("2D Neumann",   P, res, h_neu, err_neu);
        report_table("2D Robin",     P, res, h_rob, err_rob);
    }
}
