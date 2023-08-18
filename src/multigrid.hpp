#pragma once

/*
 * This is a basic geometric multigrid routine for use with LDG Poisson equation on uniform Cartesian grids.
 */

#include <unordered_map>
#include <vector>

#include "uvector.hpp"
#include "multiloop.hpp"
#include "smatrix.hpp"
#include "legendre.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include "Poisson.hpp"
#include "BlockSparseMatrix.hpp"

template<int P, int N>
class MultiGrid {

    std::vector<BlockSparseMatrix<smatrix<double, ipow(P, N)>>> Aops, Iops, Gops[N];
    std::vector<double> Mops;

    std::vector<elem_vec<P, N>> r_lev, x_lev, b_lev;

    std::vector<int> n_elements_lev;

    int levels = 1;

    int mu = 3;

public:

    MultiGrid(const uniformGrid<N> &fineGrid, const BlockSparseMatrix<smatrix<double, ipow(P, N)>> *G) {
//        levels = std::sqrt(fineGrid.get_elements_per_dim()(0)) + 1;

        int elm_dim = fineGrid.get_elements_per_dim()(0);
        while (elm_dim > 1)
        {
            elm_dim /= 2;
            ++levels;
        }

        Mops.resize(levels);
        Iops.resize(levels - 1);
        Aops.resize(levels);
        n_elements_lev.resize(levels);

        x_lev.resize(levels);
        b_lev.resize(levels);
        r_lev.resize(levels);

        for (int dim = 0; dim < N; ++dim) {
            Gops[dim].resize(levels);
        }

        n_elements_lev[0] = prod(fineGrid.get_elements_per_dim());

        for (int lev = 1; lev < levels; ++lev) {
            n_elements_lev[lev] = n_elements_lev[lev - 1] / ipow(2, N);
        }

        build_operators(fineGrid, G);

        print_operators();

    }

    int print_levels()
    {
        return levels;
    }

    elem_vec<P, N> solve(const elem_vec<P, N> &rhs) {
        b_lev[0] = rhs;
        v_cycle(0);

        std::cout << "x: " <<std::endl;
        for (int i = 0; i < n_elements_lev[0]; ++i) {
            for (int j = 0; j < ipow(P, N); ++j) {
                std::cout << x_lev[0][i](j) << std::endl;
            }
        }
        return x_lev[0];
    }

    // standard multigrid v_cycle with Gauss-Seidel smoothing
    void v_cycle(int lev) {
        if (lev < levels - 1) {
//            block_Gauss_Seidel<P, N>(Aops[lev], x_lev[lev], b_lev[lev], n_elements_lev[lev], mu);
            compute_residual(lev);
            restrict_r(lev);
            v_cycle(lev + 1);
            interpolate_x(lev);
//            block_Gauss_Seidel<P, N>(Aops[lev], x_lev[lev], b_lev[lev], n_elements_lev[lev], mu);
        } else {
            // bottom level direct solve
            smatrix<double, ipow(P, N)> Ap;
            Ap = pseudo_inverse_with_Eigen(Aops[levels - 1](0, 0));
            x_lev[levels - 1][0] = matvec(Ap, b_lev[levels - 1][0]);
        }
    }

    void print_operators() {
        std::cout << "\n--- Printing Iops ---" << std::endl;
        for (int lev = 0; lev < levels - 1; ++lev) {
            std::cout << "\n\n ---------------- Level: " << lev << " ----------------" << std::endl;
            for (int i = 0; i < n_elements_lev[lev]; ++i) {
                for (auto k : Iops[lev].row[i]) {
                    std::cout << "(" << i << "," << k << ") " << std::endl;
                    Iops[lev](i,k).print();
                }
            }
        }

        std::cout << "\n--- Printing Gops ---" << std::endl;
        for (int dim = 0; dim < N; ++dim) {
            std::cout << "Dim: " << dim << std::endl;
            for (int lev = 0; lev < levels; ++lev) {
                std::cout << "\n\n ---------------- Level: " << lev << " ----------------" << std::endl;
                for (int i = 0; i < n_elements_lev[lev]; ++i) {
                    for (auto j: Gops[dim][lev].row[i]) {
                        if (lev == levels-1){
                            std::cout << "(" << i << "," << j << ") " << std::endl;
                            Gops[dim][lev](i, j).print();
                            std::cout << std::endl;
                        }
                    }
                }
            }
        }

        std::cout << "\n--- Printing Aops ---" << std::endl;
        for (int lev = 0; lev < levels; ++lev) {
            std::cout << "Level: " << lev << std::endl;
            for (int i = 0; i < n_elements_lev[lev]; ++i) {
                for (auto j: Aops[lev].row[i]) {
                    if (lev == levels-1){
                        std::cout << "(" << i << "," << j << ") " << std::endl;
                        Aops[lev](i, j).print();
                        std::cout << std::endl;
                    }
                }
            }
        }

    }

    void build_operators(const uniformGrid<N> &fineGrid, const BlockSparseMatrix<smatrix<double, ipow(P, N)>> *G) {
        // build interpolation operator and mass matrix coefficient
        for (int lev = 0; lev < levels - 1; ++lev) {
            uniformGrid<N> grid_lev(fineGrid.get_elements_per_dim() / ipow(2, lev), fineGrid.get_xmin(),
                                    fineGrid.get_xmax());

            Mops[lev] = prod(grid_lev.get_dx());
//            std::cout << "Building Interpolation Operator on: " << lev << std::endl;
            build_interpolation_operator(grid_lev, Iops[lev]);
        }

        Mops[levels - 1] = prod(fineGrid.get_xmax() - fineGrid.get_xmin());

        // build Gops
        for (int dim = 0; dim < N; ++dim) {
            Gops[dim][0] = G[dim];
            for (int lev = 1; lev < levels; ++lev) {

                // compute (G_f * I_cf)
                BlockSparseMatrix<smatrix<double, ipow(P, N)>> GI;
                GI = BlockSparse_matmat<P, N>(Gops[dim][lev - 1], Iops[lev - 1], n_elements_lev[lev - 1]);

                // compute I_cf^T * (G_f * I_cf)
                Gops[dim][lev] = BlockSparse_matmat<P, N>(Iops[lev - 1], GI, n_elements_lev[lev], true, Mops[lev-1] / Mops[lev]);
            }
        }

        // build Aops
        for (int lev = 0; lev < levels; ++lev) {
            for (int dim = 0; dim < N; ++dim) {
                for (int i = 0; i < n_elements_lev[lev]; ++i) {
                    for (auto k: Gops[dim][lev].col[i]) {
                        for (auto j: Gops[dim][lev].row[k]) {
                            Aops[lev](i, j) +=
                                    Mops[lev] * matmat(Gops[dim][lev](k, i).transpose(), Gops[dim][lev](k, j));
                        }
                    }
                }
            }
        }
    }


    smatrix<double, ipow(P, N)> transform(const algoim::uvector<algoim::uvector<double, N>, 2> &rect_s,
                                          const algoim::uvector<algoim::uvector<double, N>, 2> &rect_d) {
        smatrix<double, ipow(P, N)> C;

        algoim::uvector<double, ipow(P, N)> basis_s, basis_d;
        algoim::uvector<double, N> pos_d, pos_s;

        double weights;

        for (algoim::MultiLoop<N> i(0, P); ~i; ++i) {
            weights = 1.;

            for (int dim = 0; dim < N; ++dim) {
                weights *= GaussQuad::w(P, i(dim));
                pos_d(dim) = GaussQuad::x(P, i(dim));
            }

            for (int dim = 0; dim < N; ++dim) {
                pos_s(dim) = (pos_d(dim) * (rect_d(1)(dim) - rect_d(0)(dim)) + rect_d(0)(dim) - rect_s(0)(dim)) /
                             (rect_s(1)(dim) - rect_s(0)(dim));
            }

            compute_basis_at_point<P, N>(basis_s, pos_s);
            compute_basis_at_point<P, N>(basis_d, pos_d);

            C += weights * outer_prod<ipow(P, N)>(basis_d, basis_s);
        }
        return C;
    }

    void build_interpolation_operator(const uniformGrid<N> &grid,
                                      BlockSparseMatrix<smatrix<double, ipow(P, N)>> &I_cf) {
        // create a coarse grid with a scaling factor of 2
        algoim::uvector<int, N> coarse_grid_elements = 0;

        for (int dim = 0; dim < N; ++dim) {
            coarse_grid_elements(dim) = int(grid.get_elements_per_dim()(dim) / 2);
        }

        uniformGrid<N> coarseGrid(coarse_grid_elements, grid.get_xmin(), grid.get_xmax());

        for (algoim::MultiLoop<N> i(0, grid.get_elements_per_dim()); ~i; ++i) {
            int source_id, dest_id;
            algoim::uvector<int, N> source_element;
            algoim::uvector<algoim::uvector<double, N>, 2> source_rect, dest_rect;

            dest_id = grid.get_element_id(i());

            for (int dim = 0; dim < N; ++dim) {
                source_element(dim) = int(i(dim) / 2);
            }

            source_id = coarseGrid.get_element_id(source_element);

            source_rect(0) = coarseGrid.get_element_min(source_element);
            source_rect(1) = source_rect(0) + coarseGrid.get_dx();

            dest_rect(0) = grid.get_element_min(i());
            dest_rect(1) = dest_rect(0) + grid.get_dx();

            I_cf(dest_id, source_id) = transform(source_rect, dest_rect);

        }
    }

    void compute_residual(const int lev) {

        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            for (auto j : Aops[lev].row[i]) {
                r_lev[lev][i] = b_lev[lev][i] - matvec(Aops[lev](i, j), x_lev[lev][j]);
            }
        }
    }

    void restrict_r(const int lev) {

        for (int i = 0; i < n_elements_lev[lev+1]; ++i) {
            for (auto j: Iops[lev].col[i]) {
                b_lev[lev + 1][i] = matvec(Iops[lev](j, i).transpose(), r_lev[lev][j]);
            }
        }
    }

    void interpolate_x(const int lev) {

        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            for (auto j: Iops[lev].row[i]) {
                x_lev[lev][i] += matvec(Iops[lev](i, j), x_lev[lev + 1][j]);
            }
        }
    }

    void print_r_lev(const int lev) {
        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            std::cout << r_lev[lev][i] << std::endl;
        }
    }

    void print_x_lev(const int lev) {
        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            std::cout << x_lev[lev][i] << std::endl;
        }
    }

    void print_b_lev(const int lev) {
        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            std::cout << b_lev[lev][i] << std::endl;
        }
    }

    void print_Aops_lev(const int lev) {
        std::cout << "\n --- Printing Aops[" << lev << "] --- " << std::endl;
        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            for (auto j: Aops[lev].row[i]) {
                std::cout << i << "," << j << std::endl;
                Aops[lev](i, j).print();
            }
        }
    }

    void print_Gops_lev(const int lev) {
        std::cout << "\n --- Printing Gops[" << lev << "] --- " << std::endl;
        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            for (auto j: Gops[0][lev].row[i]) {
                std::cout << i << "," << j << std::endl;
                Gops[0][lev](i, j).print();
            }
        }
    }


};