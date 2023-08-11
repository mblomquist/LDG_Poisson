//
// Created by mblomquist on 8/3/23.
//

#ifndef LDG_POISSON_MULTIGRID_HPP
#define LDG_POISSON_MULTIGRID_HPP

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
class MultiGrid
{

    std::vector<BlockSparseMatrix<smatrix<double, ipow(P,N)>>> Aops, Iops, Gops[N];
    std::vector<double> Mops;

    std::vector<elem_vec<P,N>> r_lev, x_lev, rhs_lev;

    std::vector<int> n_elements_lev;

    int levels = 1;



public:

    MultiGrid(const uniformGrid<N> &fineGrid, BlockSparseMatrix<smatrix<double, ipow(P, N)>> *G)
    {
        levels = std::sqrt(fineGrid.get_elements_per_dim()(0)) + 1;

        Mops.resize(levels);
        Iops.resize(levels-1);
        Aops.resize(levels);

        r_lev.resize(levels);

        for (int dim = 0; dim < N; ++dim) {
            Gops[dim].resize(levels);
        }

        for (int lev = 0; lev < levels; ++lev) {
            n_elements_lev[lev] = fineGrid.get_elements_per_dim()(0) / ipow(2, lev);
        }

        build_operators(fineGrid, G);

    }

    void build_operators(const uniformGrid<N> &fineGrid, BlockSparseMatrix<smatrix<double, ipow(P, N)>> *G)
    {
        // build interpolation operator and mass matrix coefficient
        for (int lev = 0; lev < levels-1; ++lev) {
            uniformGrid<N> grid_lev(fineGrid.get_elements_per_dim() / ipow(2, lev), fineGrid.get_xmin(), fineGrid.get_xmax());

            Mops[lev] = prod(grid_lev.get_dx());
            build_interpolation_operator(grid_lev, Iops[lev]);
        }

        Mops[levels-1] = prod(fineGrid.get_xmax() - fineGrid.get_xmin());

        // build Gops
        for (int dim = 0; dim < N; ++dim) {
            for (int lev = 1; lev < levels; ++lev) {

                // compute (G_f * I_cf)
                BlockSparseMatrix<smatrix<double, ipow(P,N)>> GI;
                for (int i = 0; i < n_elements_lev[lev]; ++i) {
                    for (auto k : Gops[dim][lev-1].row[i]){
                        for (auto j : Iops[lev-1].row[i]){
                            GI(i,j) += matmat(Gops[dim][lev-1](i,k),Iops[lev-1](k,j));
                        }
                    }
                }

                // compute I_cf^T * (G_f * I_cf)
                for (int i = 0; i < n_elements_lev[lev]; ++i) {
                    for (auto k : Gops[dim][lev-1].row[i]){
                        for (auto j : Iops[lev-1].row[i]){
                            Gops[dim][lev](k,j) += matmat(Iops[lev-1](k,i).transpose(),GI(i,j));
                        }
                    }
                }

            }
        }

        // build Aops
        for (int lev = 0; lev < levels; ++lev) {
            for (int dim = 0; dim < N; ++dim) {
                for (int i = 0; i < n_elements_lev[lev]; ++i) {
                    for (auto k : Gops[dim][lev].row[i]){
                        for (auto j : Gops[dim][lev].row[i]){
                            Aops[lev](i,j) += Mops[lev] * matmat(Gops[dim][lev](k,i).transpose(),Gops[dim][lev](k,j));
                        }
                    }
                }
            }
        }
    }


    smatrix<double, ipow(P,N)> transform(algoim::uvector<algoim::uvector<double, N>, 2> rect_s,
                                         algoim::uvector<algoim::uvector<double, N>, 2> rect_d)
    {
        GaussQuad quad;

        smatrix<double, ipow(P,N)> C;

        algoim::uvector<double, ipow(P,N)> basis_s, basis_d;
        algoim::uvector<double, N> pos_d, pos_s;

        double weights;

        for (algoim::MultiLoop<N> i(0, P); ~i; ++i)
        {
            weights = 1.;

            for (int dim = 0; dim < N; ++dim) {
                weights *= quad.w(P,i(dim));
                pos_d(dim) = quad.x(P, i(dim));
            }

            for (int dim = 0; dim < N; ++dim) {
                pos_s(dim) = (pos_d(dim) * (rect_d(1)(dim) - rect_d(0)(dim)) + rect_d(0)(dim) - rect_s(0)(dim)) / (rect_s(1)(dim) - rect_s(0)(dim));
            }

            compute_basis_at_point<P,N>(basis_s, pos_s);
            compute_basis_at_point<P,N>(basis_d, pos_d);

            C += weights * outer_prod<ipow(P,N)>(basis_d, basis_s);
        }
        return C;
    }

    void build_interpolation_operator(uniformGrid<N> &grid,
                                      BlockSparseMatrix<smatrix<double, ipow(P,N)>> &I_cf)
    {
        // create a coarse grid with a scaling factor of 2
        algoim::uvector<int, N> coarse_grid_elements = 0;

        for (int dim = 0; dim < N; ++dim) {
            coarse_grid_elements(dim) = int(grid.get_elements_per_dim()(dim) / 2);
        }

        uniformGrid<N> coarseGrid(coarse_grid_elements, grid.get_domain_min(), grid.get_domain_max());

        for (algoim::MultiLoop<N> i(0,grid.get_elements_per_dim()); ~i; ++i)
        {
            int source_id, dest_id;
            algoim::uvector<int, N> source_element;
            algoim::uvector<algoim::uvector<double, N>, 2> source_rect, dest_rect;

            dest_id = grid.get_element_id(i());

            for (int dim = 0; dim < N; ++dim) {
                source_element(dim) = int(i(dim)/2);
            }

            source_id = coarseGrid.get_element_id(source_element);

            source_rect(0) = coarseGrid.get_element_min(source_element);
            source_rect(1) = source_rect(0) + coarseGrid.get_dx();

            dest_rect(0) = grid.get_element_min(i());
            dest_rect(1) = dest_rect(0) + grid.get_dx();

            I_cf(dest_id, source_id) = transform(source_rect, dest_rect);
        }
    }

    void compute_residual(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &x_lev,
                          std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &f_lev,
                          int lev)
    {
        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            for (auto j : Aops[lev].row[i]) {
                r_lev[lev][i] += f_lev - matvec(Aops[lev](i,j),x_lev[j]);
            }
        }
    }

    void restrict_r(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &f_lev,
                    int lev)
    {

        for (int i = 0; i < n_elements_lev[lev]; ++i) {
            for (auto j : Aops[lev].row[i]) {
                r_lev[lev+1][j] += matvec(Iops[lev](j,i).transpose(),f_lev[i]);
            }
        }
    }
    
    void vcycle(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &x_lev,
                std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &f_lev,
                int lev)
    {
        if (lev < levels-1)
        {
            block_Gauss_Seidel<P,N>(Aops[lev], x_lev, f_lev, n_elements_lev[lev]);
            compute_residual(x_lev, f_lev);
            restrict_r(r_lev[lev], lev);
        }
        else
        {
            // direct solve

            smatrix<double, ipow(P,N)> Ap;
            Ap = pseudo_inverse_with_Eigen(Aops[levels-1]);
            x_lev[0] = matvec(Ap,f_lev[0]);
        }
    }
};





#endif //LDG_POISSON_MULTIGRID_HPP
