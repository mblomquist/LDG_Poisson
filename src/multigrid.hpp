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

    std::vector<BlockSparseMatrix<smatrix<double, ipow(P,N)>>> Aops, Gops, Iops;
    std::vector<double> Mops;

    int levels = 1;

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
//            pos_s(dim) = (pos_d(dim) * (rect_s(1)(dim) - rect_s(0)(dim)) + rect_s(0)(dim) - rect_d(0)(dim)) / (rect_d(1)(dim) - rect_d(0)(dim));
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

public:

    MultiGrid(uniformGrid<N> &fineGrid, BlockSparseMatrix<smatrix<double, ipow(P, N)>> G[N])
    {
        levels = fineGrid.get_elements_per_dim()(0) / 2 + 1;

        for (int lev = 0; lev < levels; ++lev) {
            uniformGrid<N> grid_lev(fineGrid.get_elements_per_dim() / ipow(2, lev), fineGrid.get_xmin(), fineGrid.get_xmax());

            Mops[lev] = prod(grid_lev.get_dx());
            build_interpolation_operator(grid_lev, Iops[lev]);

        }
    }
};


template<int P, int N>
void block_Gauss_Seidel(BlockSparseMatrix<smatrix<double, ipow(P, N)>> &A,
                        std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &x,
                        std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &b,
                        int num_elements,
                        int max_itr = 3)
{
    for (int itr = 0; itr < max_itr; ++itr)
    {
        for (int i = 0; i < num_elements; ++i)
        {
            algoim::uvector<double, ipow(P,N)> sum = 0.;

            for (auto j : A.row[i])
            {
                if (i != j)
                    sum += matvec(A(i,j),x[j]);
            }

            x[i] = matvec(pseudo_inverse_with_Eigen(A(i, i)), b[i] - sum);
        }
    }
}


#endif //LDG_POISSON_MULTIGRID_HPP
