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

    L_f_ij = (c1-1.)*A_ii + c2*A_ij - c1*A_ij - (c2-1.)*A_jj;

    std::cout << "\nL_f_ij:" << std::endl;
    print_small_matrix<P,N>(L_f_ij);

}


template<int P, int N>
algoim::uvector<smatrix<double, ipow(P,N)>, 4> compute_lifting_operator_on_ref_face(int dim)
{
    GaussQuad quad;
    constexpr int Q = int((2*P+1)/2)+1;

    algoim::uvector<smatrix<double, ipow(P,N)>, 4> L_f_ij;

    smatrix<double, ipow(P,N)> A_ii, A_ij, A_ji, A_jj;

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

    L_f_ij(0) = A_ii;
    L_f_ij(1) = A_ij;
    L_f_ij(2) = A_ji;
    L_f_ij(3) = A_jj;

    return L_f_ij;
}

struct HashTuple{
    std::size_t operator()(const std::tuple<int, int>& t) const {
        std::size_t hash_0 = std::hash<int>()(std::get<0>(t));
        std::size_t hash_1 = std::hash<int>()(std::get<1>(t));
        return hash_0 ^ (hash_1 << 1);
    }
};

struct KeyTupleEqual{
    bool operator()(const std::tuple<int, int>& left_t, const std::tuple<int, int>& right_t) const {
        return std::equal_to<int>()(std::get<0>(left_t), std::get<0>(right_t)) &&
               std::equal_to<int>()(std::get<1>(left_t), std::get<1>(right_t));
    }
};

template<int P, int N>
void compute_lifting_operator_periodic_grid(uniformGrid<N> grid)
{

    std::unordered_map<std::tuple<int,int>, smatrix<double, ipow(P,N)>, HashTuple, KeyTupleEqual> L[N];
    algoim::uvector<smatrix<double, ipow(P,N)>, 4> L_f_ij;

    algoim::uvector<int, N> elements_per_dim = grid.get_elements_per_dim();

    double c1 = 1.;
    double c2 = 1.-c1;

    for (int dim = 0; dim < N; ++dim) {

        for (MultiLoop<N> i(0,elements_per_dim); ~i; ++i)
        {
            algoim::uvector<int, N> element_i, element_j;

            element_i(dim) = (i(dim) - 1 == -1) ? elements_per_dim(dim) - 1 : i(dim) - 1;
            element_j(dim) = i(dim);

            for (int d = 0; d < N; ++d) {
                if (d != dim){
                    element_i(d) = i(d);
                    element_j(d) = i(d);
                }
            }

            L_f_ij = compute_lifting_operator_on_ref_face<P,N>(dim);

            L[dim][{grid.get_element_id(element_i), grid.get_element_id(element_i)}] = (c1-1.)*L_f_ij(0);
            L[dim][{grid.get_element_id(element_i), grid.get_element_id(element_j)}] = c2*L_f_ij(1);
            L[dim][{grid.get_element_id(element_j), grid.get_element_id(element_i)}] = -c1*L_f_ij(2);
            L[dim][{grid.get_element_id(element_j), grid.get_element_id(element_j)}] = (1.-c2)*L_f_ij(3);

        }
    }
}

#endif //LDG_POISSON_SCRATCH_HPP
