//
// Created by mblomquist on 7/25/23.
//

#ifndef LDG_POISSON_SCRATCH_HPP
#define LDG_POISSON_SCRATCH_HPP

#include <iostream>
#include <functional>
#include "quadrature.hpp"
#include "uvector.hpp"
#include "grid.hpp"
#include "Poisson.hpp"

using namespace algoim;

void suppressed_multiloop()
{

    GaussQuad quad;

    constexpr int Q = 3; // quad points
    constexpr int D = 3; // dimension

    uvector<int, D> min, max;

    min = 0;
    max = Q;

    uvector<double, D> pos;
    uvector<double, D-1> pos_Dmo;

    for (int dim = 0; dim < D; ++dim) {

        std::cout << "\nFace: " << dim << std::endl;
        pos = 0;
        pos(dim) = 1.;

        for (MultiLoop<D-1> i(0,Q); ~i; ++i)
        {
            for (int j = 0; j < D-1; ++j) {
                    pos_Dmo(j) = quad.x(Q,i(j));
            }

            int t = 0;
            for (int j = 0; j < D; ++j) {
                if (j != dim){
                    pos(j) = pos_Dmo(t);
                    ++t;
                }
            }
            std::cout << dim << " " << i() << " " << pos << std::endl;
        }
    }
}

template<int N>
double evaluate_face_integral(std::function<double(algoim::uvector<double, N> x)> u,
                              std::function<double(algoim::uvector<double, N> x)> v)
{
    GaussQuad quad;

    constexpr int Q = 3;

    uvector<double, N> eval_pos_u, eval_pos_v;
    uvector<double, N-1> pos_Dmo;

    eval_pos_u(0) = 1.; // u is on the left
    eval_pos_v(0) = 0.; // v is on the right

    int dim = 0;

    double sum = 0.;

    for (MultiLoop<N-1> i(0,Q); ~i; ++i)
    {
        double weight = 1.;

        for (int j = 0; j < N-1; ++j) {
            pos_Dmo(j) = quad.x(Q,i(j));
            weight *= quad.w(Q,i(j));
        }

        int t = 0;
        for (int j = 0; j < N; ++j) {
            if (j != dim){
                eval_pos_u(j) = pos_Dmo(j);
                eval_pos_v(j) = pos_Dmo(j);
                ++t;
            }

        }

        sum += weight * u(eval_pos_u) * v(eval_pos_v);
    }

    return sum;
}

template<int P, int N>
void print_smatrix(smatrix<double, ipow(P,N)> A)
{
    for (int i = 0; i < ipow(P, N); ++i) {
        for (int j = 0; j < ipow(P, N); ++j) {
            std::cout << A(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template<int P, int N>
void compute_basis_quadrature()
{
    GaussQuad quad;
    constexpr int Q = int((2.*P+1.)/2.)+1;

    smatrix<double, ipow(P,N)> A_ii, A_ij, A_ji, A_jj;

    uvector<double, ipow(P,N)> eval_i, eval_j;

    uvector<double, N> eval_pos_i, eval_pos_j;
    uvector<double, N-1> pos_Dmo;

    int dim = 0;

    eval_pos_i(dim) = 1.; // element to the left of the face
    eval_pos_j(dim) = 0.; // element to the right of the face

    for (MultiLoop<N - 1> i(0, Q); ~i; ++i) {
        double weight = 1.;

        for (int j = 0; j < N - 1; ++j) {
            pos_Dmo(j) = quad.x(Q, i(j));
            weight *= quad.w(Q, i(j));
        }

        int t = 0;
        for (int j = 0; j < N; ++j) {
            if (j != dim) {
                eval_pos_i(j) = pos_Dmo(j);
                eval_pos_j(j) = pos_Dmo(j);
                ++t;
            }

        }

        eval_i = evaluate_basis_at_point_on_face<P, N>(eval_pos_i);
        eval_j = evaluate_basis_at_point_on_face<P, N>(eval_pos_j);

        A_ii += outer_prod(eval_i, eval_i);
        A_ij += outer_prod(eval_i, eval_j);
        A_ji += outer_prod(eval_j, eval_i);
        A_jj += outer_prod(eval_j, eval_j);
    }

    std::cout << "\nA_ii:" << std::endl;
    print_smatrix<P,N>(A_ii);

    std::cout << "\nA_ij:" << std::endl;
    print_smatrix<P,N>(A_ij);

    std::cout << "\nA_ji:" << std::endl;
    print_smatrix<P,N>(A_ji);

    std::cout << "\nA_jj:" << std::endl;
    print_smatrix<P,N>(A_jj);

}



#endif //LDG_POISSON_SCRATCH_HPP
