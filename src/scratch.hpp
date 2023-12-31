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
void print_small_matrix(smatrix<double, ipow(P,N)> A)
{
    for (int i = 0; i < ipow(P, N); ++i) {
        for (int j = 0; j < ipow(P, N); ++j) {
            if (std::abs(A(i,j)) > 1.0e-12)
                std::cout << A(i, j) << " ";
            else
                std::cout << "0 ";
        }
        std::cout << std::endl;
    }
}

template<int P, int N>
void compute_lifting_operator_on_a_face(int dim)
{
    GaussQuad quad;
    constexpr int Q = int((2*P+1)/2)+1;

    smatrix<double, ipow(P,N)> A_ii, A_ij, A_ji, A_jj, L_f_ij;

    uvector<double, ipow(P,N)> eval_i, eval_j;

    uvector<double, N> eval_pos_i, eval_pos_j;
    uvector<double, N-1> pos_Dmo;

    double c1 = 1.;
    double c2 = 1.-c1;

    eval_pos_i(dim) = 1.;
    eval_pos_j(dim) = 0.;

    for (MultiLoop<N-1> i(0, Q); ~i; ++i) {
        double weight = 1.;

        for (int j = 0; j < N - 1; ++j) {
            pos_Dmo(j) = quad.x(Q, i(j));
            weight *= quad.w(Q, i(j));
        }

        int t = 0;
        for (int j = 0; j < N; ++j) {
            if (j != dim) {
                eval_pos_i(j) = pos_Dmo(t);
                eval_pos_j(j) = pos_Dmo(t);
                ++t;
            }

        }

        eval_i = evaluate_basis_as_point<P, N>(eval_pos_i);
        eval_j = evaluate_basis_as_point<P, N>(eval_pos_j);

        A_ii += outer_prod(eval_i, eval_i) * weight;
        A_ij += outer_prod(eval_i, eval_j) * weight;
        A_ji += outer_prod(eval_j, eval_i) * weight;
        A_jj += outer_prod(eval_j, eval_j) * weight;

    }

}

template<int P, int N>
void test_transform()
{
    algoim::uvector<algoim::uvector<double, N>, 2> source, dest;
    smatrix<double, ipow(P,N)> C_sd;

    source(0) = 10.0;
    source(1) = 11.0;

    dest(0) = 10.5;
    dest(1) = 12.0;

    C_sd = transform<P,N>(source, dest);

    constexpr int num_pts = 5;

    algoim::uvector<algoim::uvector<double, N>, num_pts> rnd_pts, x_s, x_d;

    for (int i = 0; i < num_pts; ++i) {
        for (int dim = 0; dim < N; ++dim) {
            rnd_pts(i)(dim) = static_cast<double>(std::rand()) / RAND_MAX;
        }
    }

    for (int i = 0; i < num_pts; ++i) {
        for (int dim = 0; dim < N; ++dim) {
            x_s(i)(dim) = (rnd_pts(i)(dim) - source(0)(dim))/ (source(1)(dim) - source(0)(dim));
            x_d(i)(dim) = (rnd_pts(i)(dim) - dest(0)(dim)) / (dest(1)(dim)   - dest(0)(dim));
        }
    }

    algoim::uvector<double, ipow(P,N)> u, Cu, Ctu;
    algoim::uvector<double, num_pts> result0, result1;

    for (int i = 0; i < ipow(P,N); ++i) {
        u(i) = static_cast<double>(std::rand()) / RAND_MAX;
    }

    Cu = matvec<ipow(P,N)>(C_sd, u);

    std::cout << "\n--- Random Points ---" << std::endl;
    for (int i = 0; i < num_pts; ++i) {
        std::cout << "pt " << i << ": " << rnd_pts(i) << ", " << x_s(i) << ", " << x_d(i) << std::endl;
    }

    std::cout << "\n--- Results ---" << std::endl;
    std::cout << "source rectangle: " << source << std::endl;
    std::cout << "destination rectangle: " << dest << std::endl;
    for (int i = 0; i < num_pts; ++i) {
        result0(i) = compute_basis_coefficients_at_point<P,N>(u,x_s(i));
        result1(i) = compute_basis_coefficients_at_point<P,N>(Cu,x_d(i));

        std::cout << "pt " << i << ": " << std::abs(result0(i) - result1(i)) << std::endl;
    }
}



#endif //LDG_POISSON_SCRATCH_HPP
