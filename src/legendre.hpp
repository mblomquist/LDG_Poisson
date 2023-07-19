//
// Created by mblomquist on 7/7/23.
//

#ifndef DG_UTILS_LEGENDRE_HPP
#define DG_UTILS_LEGENDRE_HPP

#include "uvector.hpp"
#include "multiloop.hpp"
#include "quadrature.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include <cmath>
#include <vector>
#include <functional>

template<int P, int N>
int unfold(const algoim::uvector<int, N> indx_t)
{
    int indx_f = 0;

    for (int i = 0; i < N; ++i)
    {
        indx_f += indx_t(i) * ipow(P,i);
    }

    return indx_f;
}

template<int P>
void evaluate_legendre_at_point(const double x, algoim::uvector<double, P> &px)
{
    px(0) = 1.;

    if constexpr (P >= 1)
        px(1) = x;

    for (int i = 1; i < P-1; ++i) {
        px(i+1) = ((2.*i+1.)*x*px(i) - i*px(i-1))/(i+1.);
    }
}

template<int P>
void evaluate_shifted_legendre(const double x, algoim::uvector<double, P> &px)
{
    evaluate_legendre_at_point<P>((2.*x-1.), px);

    // normalize
    for (int i = 0; i < P; ++i) {
        px(i) *= std::sqrt(2.*i+1.);
    }
}

template<int P, int N>
void compute_basis_at_point(algoim::uvector<double, ipow(P,N)> &basis_Nd,
                            algoim::uvector<double, N> &pt)
{
    algoim::uvector<algoim::uvector<double,P>, N> basis_1d = 0.;

    for (int dim = 0; dim < N; ++dim) {
        evaluate_shifted_legendre<P>(pt(dim), basis_1d(dim));
    }

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
    {
        double prod = 1.;
        for (int dim = 0; dim < N; ++dim) {
            prod *= basis_1d(dim)(i(dim));
        }
        basis_Nd(unfold<P,N>(i())) = prod;
    }

}

template<int P, int N>
void l2_projection_on_reference_element(const std::function<double(const algoim::uvector<double,N>&)> func,
                                 algoim::uvector<double, ipow(P,N)> &p_fun)
{
    GaussQuad quad;
    constexpr int Q = int((2.*P+1.)/2.)+1;

    // create a stacked variable of the basis at multiple quadrature points
    algoim::uvector<algoim::uvector<double, ipow(P,N)>, ipow(Q,N)> fat_basis;

    algoim::uvector<double, N> pos;
    double weights;

    for (algoim::MultiLoop<N> i(0, Q); ~i; ++i)
    {
        weights = 1.;

        for (int dim = 0; dim < N; ++dim) {
            weights *= quad.w(Q,i(dim));
            pos(dim) = quad.x(Q,i(dim));
        }

        compute_basis_at_point<P,N>(fat_basis(unfold<Q,N>(i())), pos);

        p_fun += fat_basis(unfold<Q,N>(i())) * weights * func(pos);

    }

}

template<int P, int N>
double evaluate_basis_at_point(algoim::uvector<double, ipow(P,N)> coeff,
                               algoim::uvector<double, N> pt)
{
    algoim::uvector<algoim::uvector<double,P>, N> basis = 0.;

    for (int dim = 0; dim < N; ++dim) {
        evaluate_shifted_legendre<P>(pt(dim), basis(dim));
    }

    double sum = 0.;

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
    {
        double prod = 1.;
        for (int dim = 0; dim < N; ++dim) {
            prod *= basis(dim)(i(dim));
        }
        sum += coeff(unfold<P,N>(i())) * prod;
    }

    return sum;
}


#endif //DG_UTILS_LEGENDRE_HPP